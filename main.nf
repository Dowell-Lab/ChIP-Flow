#!/usr/bin/env nextflow
/*
========================================================================================
                         ChIP-Flow - ChIP-seq PIPELINE
========================================================================================
 ChIP-seq Analysis Pipeline. Started 2018-11-06.
 #### Homepage / Documentation
https://github.com/Dowell-Lab/ChIP-Flow
 #### Authors
 Margaret Gruca <magr0763@colorado.edu>, Lynn Sanford <lysa8537@colorado.edu>
========================================================================================
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq
        SRA tools -- fastq-dump sra to generate fastq file
        FastQC (pre-trim) -- perform pre-trim FastQC on fastq files

    2. Trimming & Mapping
        BBDuk -- trim fastq files for quality and adapters
        FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)
        Mapping w/ HISAT2 -- map to genome reference file

    3. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM, BAM --> CRAM, index CRAM

    4. Quality control
        preseq -- estimate library complexity
        RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
        Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/neg reads, intron/exon ratio
        Picard -- Mark duplicates, GC content, deduplicate

    5. Visualization
        BEDTools -- Generate bedGraphs
        kentUtils -- Generate bigWigs
        IGV Tools -- Generate TDF from bedGraph

    6. MultiQC -- Generate QC report for pipeline

    7. Pipeline report


*/


def helpMessage() {
    log.info"""
    =========================================
     ChIP-Flow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile fiji --fastqs '/project/*_{R1,R2}*.fastq' --outdir '/project/'

    Required arguments:
         -profile                      Configuration profile to use. <base, fiji>
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.

    Performance options:
        --threadfqdump                 Runs multi-threading for fastq-dump for sra processing.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.
        --saveAll                      Compresses and saves all fastq reads.
        --skipBAM                      Skip saving BAM files.
        --savedup                      Save deduplicated/marked duplicate BAM files (using picard, cannot be used with --skippicard).

    QC Options:
        --skipMultiQC                  Skip running MultiQC.
        --skipFastQC                   Skip running FastQC.
        --skipRSeQC                    Skip running RSeQC.
        --skippicard                   Skip running picard.        
        --skippreseq                   Skip running preseq.
        --skippileup                   Skip running pileup.sh.
        --skipAllQC                    Skip running all QC.
        --noTrim                       Skip trimming and map only.
        --dedup                        Remove sequencing duplicates from BAM files (using picard, cannot be used with --skippicard).
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.bbmap_adapters = "$baseDir/bin/adapters.fa"
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.rcc = "$baseDir/bin/rcc.py"
params.workdir = "./nextflowTemp"
params.extract_fastqc_stats = "$baseDir/bin/extract_fastqc_stats.sh"

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

import java.text.SimpleDateFormat
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd")
output_date =  sdf.format(date)

// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.chrom_sizes ){
    chrom_sizes = file(params.chrom_sizes)
    if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }

if ( params.bbmap_adapters ){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

if ( params.picard_path ){
    picard_path = file(params.picard_path)
}

if ( params.hisat2_indices ){
    hisat2_indices = file("${params.hisat2_indices}")
}

if ( params.genome_refseq ){
    genome_refseq = file("${params.genome_refseq}")
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



/*
 * Create a channel for input read files
 */
if (params.fastqs) {
    if (params.singleEnd) {
        fastq_reads_qc = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_hisat2_notrim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
    } else {
        Channel
            .fromFilePairs( params.fastqs, size: 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_hisat2_notrim }
    }
}

else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_hisat2_notrim }
}

if (params.sras) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
}

else {
    read_files_sra = Channel.empty()
}

// Header log info
log.info """=======================================================
ChIP-Flow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'ChIPFlow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Fastqs']           = params.fastqs
summary['SRAs']             = params.sras
summary['Genome Ref']       = params.genome
summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Save Dup']         = params.savedup ? 'YES' : 'NO'
summary['Remove Dup']       = params.dedup ? 'YES' : 'NO'
summary['Run FastQC']       = params.skipFastQC ? 'NO' : 'YES'
summary['Run preseq']       = params.skippreseq ? 'NO' : 'YES'
summary['Run pileup']       = params.skippileup ? 'NO' : 'YES'
summary['Run RSeQC']        = params.skipRSeQC ? 'NO' : 'YES'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Skip All QC']      = params.skipAllQC ? 'YES' : 'NO'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    time '1h'

    output:
    stdout into software_versions

    script:
    """
    printf "chipflow_version: %s\n" ${params.version}
    printf "nextflow_version: %s\n" ${workflow.nextflow.version}
    printf "fastqc_version: %s\n" \$(fastqc -v | head -1 | awk -F " v" '{print \$2}')
    printf "bbmap_version: %s\n" \$(bbversion.sh --version | head -1)
    printf "hisat2_version: %s\n" \$(hisat2 --version | head -1 | awk '{print \$NF}')
    printf "samtools_version: %s\n" \$(samtools --version | head -1 | awk '{print \$NF}')
    printf "picard_gc_version: %s\n" \$(java -jar ${params.picard_path} CollectGcBiasMetrics --version 2>&1 | head -1 | awk -F "-" '{print \$1}')
    printf "picard_dups_version: %s\n" \$(java -jar ${params.picard_path} MarkDuplicates --version 2>&1 | head -1 | awk -F "-" '{print \$1}')
    printf "preseq_version: %s\n" \$(preseq 2>&1 | head -2 | tail -1 | awk '{print \$NF}')
    printf "bedtools_version: %s\n" \$(bedtools --version | head -1 | awk -F " v" '{print \$2}')
    printf "igvtools_version: %s\n" \$(igvtools version | head -1 | awk '{print \$3}')
    printf "rseqc_version: %s\n" \$(infer_experiment.py --version | head -1 | awk '{print \$NF}')
    printf "pipeline_hash: %s\n" ${workflow.scriptId}
    """
}

software_versions.collectFile(name: "software_versions_chipflow_${output_date}_${workflow.runName}.yaml", storeDir: "${params.outdir}/pipeline_info")

/*
 * Step 1 -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"
    if (params.threadfqdump) {
        cpus 8 }
    else {
        cpus 1
    }
    if (params.savefq || params.saveAllfq) {
        publishDir "${params.outdir}/fastq", mode: 'copy'
    }
    
    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    tuple val(prefix), file("*.fastq.gz") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_hisat2_notrim_sra
   

    script:
    prefix = reads.baseName
    if (!params.threadfqdump && !params.singleEnd) {
        """
        echo ${prefix}
        fastq-dump --split-3 ${reads} --gzip
        """
    } else if (!params.singleEnd) {
         """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --split-3 \
            --sra-id ${reads}
        """
    } else if (!params.threadfqdump) {
        """
        echo ${prefix}
        fastq-dump ${reads} --gzip
        """
    } else {
        """
        export PATH=~/.local/bin:$PATH
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --sra-id ${reads}
        """
    }
}

/*
 * STEP 1+ - FastQC
 */

process fastQC {
    tag "$prefix"
    memory '8 GB'
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("zip") > 0)     "qc/fastqc/zips/$filename"
        else if (filename.indexOf("html") > 0)    "qc/fastqc/$filename"
        else if (filename.indexOf("txt") > 0)     "qc/fastqc_stats/$filename"
        else null            
    }
    
    when:
    !params.skipFastQC && !params.skipAllQC

    input:
    tuple val(prefix), path(reads) from fastq_reads_qc.mix(fastq_reads_qc_sra)

    output:
    file "*.{zip,html}" into fastqc_results
    file "*.fastqc_stats.txt" into fastqc_stats

    script:
    def reads_1 = ''
    def reads_2 = ''
    def prefix_1 = ''
    def prefix_2 = ''
    if (params.singleEnd) {
        reads_1 = reads.first()
        prefix_1 = reads_1.simpleName
    } else {
        reads_1 = reads.first()
        reads_2 = reads.last()
        prefix_1 = reads_1.simpleName
        prefix_2 = reads_2.simpleName
    }
    if (params.singleEnd) {
    """
    echo ${prefix_1}
    fastqc ${reads_1}
    
    ${params.extract_fastqc_stats} \
        --srr=${prefix_1} \
        > ${prefix_1}.fastqc_stats.txt
    """
    } else {
    """
    echo ${prefix_1} ${prefix_2}
    fastqc ${reads_1}

    ${params.extract_fastqc_stats} \
        --srr=${prefix_1} \
        > ${prefix_1}.fastqc_stats.txt

    fastqc ${reads_2}

    ${params.extract_fastqc_stats} \
        --srr=${prefix_2} \
        > ${prefix_2}.fastqc_stats.txt
    """
    }
}

/*
 * STEP 2 - Trimming & Mapping
 */

process bbduk_hisat2 {
    tag "$name"
    cpus 32
    memory '40 GB'
    publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*.{refstats,trimstats}.txt"
    publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*hisat2_mapstats.txt"    
    if (params.saveTrim || params.saveAllfq) {
        publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    }
    
    when:
    !params.noTrim

    input:
    file(indices) from hisat2_indices
    val(indices_path) from hisat2_indices        
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    set val(name), file("*.trim.fastq.gz") into trimmed_reads_fastqc
    file "*.{refstats,trimstats}.txt" into trim_stats
    set val(name), file("*.sam") into hisat2_sam
    file("*hisat2_mapstats.txt") into hisat2_mapstats

    script:      
    prefix_pe = reads[0].toString() - ~/(_1\.)?(_R1)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
    prefix_se = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
    
    if (!params.singleEnd) {
        """
        echo ${prefix_pe}

        bbduk.sh -Xmx40g \
                  t=32 \
                  in=${reads[0]} \
                  in2=${reads[1]} \
                  out=${prefix_pe}_1.trim.fastq.gz \
                  out2=${prefix_pe}_2.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=r trimq=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=20 \
                  tpe tbo \
                  stats=${prefix_pe}.trimstats.txt \
                  refstats=${prefix_pe}.refstats.txt
                  
        hisat2 -p 32 \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${indices_path} \
               -1 ${prefix_pe}_1.trim.fastq.gz \
               -2 ${prefix_pe}_2.trim.fastq.gz \
               --new-summary \
               > ${prefix_pe}.sam \
               2> ${prefix_pe}.hisat2_mapstats.txt                     
        """
    } else {
        """
        module load bbmap/38.05
        echo ${prefix_se}
        
        bbduk.sh -Xmx40g \
                  t=32 \
                  in=${reads} \
                  out=${prefix_se}.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=r trimq=10 k=23 mink=11 hdist=1 \
                  maq=10 minlen=20 \
                  stats=${prefix_se}.trimstats.txt \
                  refstats=${prefix_se}.refstats.txt
                  
        hisat2 -p 32 \
               --very-sensitive \
               --no-spliced-alignment \
               -x ${indices_path}\
               -U ${prefix_se}.trim.fastq.gz \
               --new-summary \
               > ${prefix_se}.sam \
               2> ${prefix_se}.hisat2_mapstats.txt                  
        """
    }
}


/*
 * STEP 2+ - Trimmed FastQC
 */

process fastqc_trimmed {
    tag "$prefix"
    memory '4 GB'
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
    
    when:
    !params.skipFastQC

    input:
    set val(prefix), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html,txt}" into trimmed_fastqc_results

    script:
    """
    echo ${prefix}

    fastqc $trimmed_reads
    """
}


/*
 * STEP 2+ - Mapping w/o trimming
 */

if (params.noTrim) {
    process hisat2 {
        tag "$name"
        cpus 32
        memory '40 GB'
        time '2h'
        publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*.txt"
    
        input:
        file(indices) from hisat2_indices
        val(indices_path) from hisat2_indices
        set val(name), file(reads) from fastq_reads_hisat2_notrim_sra.mix(fastq_reads_hisat2_notrim)
    
        output:
        set val(name), file("*.sam") into hisat2_sam
        file("*.txt") into hisat2_mapstats
    
        script:
        prefix_pe = trimmed_reads[0].toString() - ~/(_1\.)?(_R1)?(flip)?(trim)?(\.flip)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
        prefix_se = trimmed_reads[0].toString() - ~/(\.flip)?(\.fq)?(\.fastq)?(\.gz)?$/
        
        def rnastrandness = ''
        if (params.forwardStranded && !params.unStranded){
            rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
        } else if (params.reverseStranded && !params.unStranded){
            rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
        }
        
        if (!params.singleEnd) {
            """
            echo ${prefix_pe}
        
            hisat2 -p 32 \
                   --very-sensitive \
                   --no-spliced-alignment \
                    $rnastrandness \
                   -x ${indices_path} \
                   -1 ${reads[0]} \
                   -2 ${reads[1]} \
                   --new-summary \
                   > ${prefix_pe}.sam \
                   2> ${prefix_pe}.hisat2_mapstats.txt                
            """
        }
        else {
            """
            echo ${prefix_se}
        
            hisat2  -p 32 \
                    --very-sensitive \
                    --no-spliced-alignment \
                    $rnastrandness \
                    -x ${indices_path} \
                    -U ${reads} \
                    --new-summary \
                    > ${prefix_se}.sam \
                    2> ${prefix_se}.hisat2_mapstats.txt                
            """
        }
    }
}


/*
 * STEP 3 - Convert to BAM format and sort
 */

process samtools {
    tag "$name"
    memory '100 GB'
    cpus 16
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if ((filename.indexOf("sorted.bam") > 0) & !params.skipBAM)                                                                                                                             "mapped/bams/$filename"
        else if ((filename.indexOf("sorted.bam.bai") > 0) & !params.skipBAM)                                                                                                                         "mapped/bams/$filename"
        else if (filename.indexOf("flagstat") > 0)                    "qc/mapstats/$filename"
        else if (filename.indexOf("millionsmapped") > 0)              "qc/mapstats/$filename"
        else if (filename.indexOf("sorted.cram") > 0)                 "mapped/crams/$filename"
        else if (filename.indexOf("sorted.cram.crai") > 0)            "mapped/crams/$filename"
        else null            
    }

    input:
    set val(name), file(mapped_sam) from hisat2_sam

    output:
    set val(name), file("${name}.sorted.bam") into sorted_bam_ch
    set val(name), file("${name}.sorted.bam.bai") into sorted_bam_indices_ch
    set val(name), file("${name}.flagstat") into bam_flagstat
    set val(name), file("${name}.millionsmapped") into bam_milmapped_bedgraph
    set val(name), file("${name}.sorted.cram") into cram_out
    set val(name), file("${name}.sorted.cram.crai") into cram_index_out

    script:
    """

    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x904 -c ${name}.sorted.bam > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
}

sorted_bam_ch
   .into { sorted_bams_for_preseq; sorted_bams_for_rseqc; sorted_bams_for_dreg_prep; sorted_bams_for_pileup; sorted_bams_for_picard; sorted_bam_for_bedgraph }

sorted_bam_indices_ch
    .into { sorted_bam_indicies_for_pileup; sorted_bam_indices_for_preseq; sorted_bam_indices_for_rseqc; sorted_bam_indices_for_picard; sorted_bam_index_for_bedgraph }

/*
 *STEP 4+ - Picard tools
 */

process picard {
    tag "$name"
    memory '20 GB'
    time '8h'
    errorStrategy 'ignore'
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("marked_dup_metrics.txt") > 0)                                   "qc/picard/dups/$filename"
        else if (filename.indexOf("gc_bias_metrics.pdf") > 0)                                      "qc/picard/gc_bias/$filename"
        else if (filename.indexOf("gc_bias_metrics.txt") > 0)                                      "qc/picard/gc_bias/$filename"
        else if (filename.indexOf("summary_metrics.txt") > 0)                                      "qc/picard/gc_bias/$filename"
        else if ((filename.indexOf("bam") > 0) && (params.savedup))                                "mapped/picard_dedup/$filename"
        else if ((filename.indexOf("bam.bai") > 0) && (params.savedup))                            "mapped/picard_dedup/$filename"            
        else null            
    }
    
    when:
    !params.skippicard 

    input:
    set val(name), file(bam_file) from sorted_bams_for_picard
    file(bam_indices) from sorted_bam_indices_for_picard

    output:
    file "*.{txt,pdf}" into picard_stats_multiqc, picard_stats_nqc
    file "*.bam" into picard_dup_bam
    
    script:
    if (!params.dedup) {
    """
    java -jar -Xmx20g ${params.picard_path} MarkDuplicates \
         I=${bam_file} \
         O=${name}.marked_duplicates.bam \
         M=${name}.marked_dup_metrics.txt             
         
    java -jar -Xmx20g ${params.picard_path} CollectGcBiasMetrics \
          I=${bam_file} \
          O=${name}.gc_bias_metrics.txt \
          CHART=${name}.gc_bias_metrics.pdf \
          S=${name}.summary_metrics.txt \
          R=${genome}    
    """
    } else {
    """
    java -jar -Xmx20g ${params.picard_path} MarkDuplicates \
         I=${bam_file} \
         REMOVE_SEQUENCING_DUPLICATES=true \
         O=${name}.dedup_unsorted.bam \
         M=${name}.marked_dup_metrics.txt
    samtools sort -@ 16 ${name}.dedup_unsorted.bam > ${name}.dedup.bam
    rm ${name}.dedup_unsorted.bam
    samtools index ${name}.dedup.bam ${name}.dedup.bam.bai
         
    java -jar -Xmx20g ${params.picard_path} CollectGcBiasMetrics \
          I=${bam_file} \
          O=${name}.gc_bias_metrics.txt \
          CHART=${name}.gc_bias_metrics.pdf \
          S=${name}.summary_metrics.txt \
          R=${genome}    
    """        
    }
}
/*
 *STEP 4+ - Preseq
 */

process preseq {
    tag "$name"
    memory '20 GB'
    time '8h'
    publishDir "${params.outdir}/qc/preseq/", mode: 'copy', pattern: "*.txt"
    
    when:
    !params.skippreseq

    input:
    set val(name), file(bam_file) from sorted_bams_for_preseq
    file(bam_indices) from sorted_bam_indices_for_preseq

    output:
    file("*.txt") into preseq_results

    script:
    """
    preseq c_curve -B -o ${name}.c_curve.txt \
           ${bam_file}

    preseq lc_extrap -B -o ${name}.lc_extrap.txt \
           ${bam_file}
    """
 }


/*
 *STEP 4+ - Analyze read distributions using RSeQC
 */

process rseqc {
    tag "$name"
    time '8h'
    memory '40 GB'
    publishDir "${params.outdir}/qc/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else null
        }

    input:
    set val(name), file(bam_file) from sorted_bams_for_rseqc
    file(bam_indices) from sorted_bam_indices_for_rseqc

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """
    read_distribution.py -i ${bam_file} \
                         -r ${genome_refseq} \
                         > ${name}.read_dist.txt

    read_duplication.py -i ${bam_file} \
                        -o ${name}.read_duplication

    infer_experiment.py -i ${bam_file} \
                        -r ${genome_refseq} \
                        > ${name}.infer_experiment.txt
    """
 }



/*
 *STEP 4+ - Analyze coverage using pileup.sh
 */

process pileup {
    tag "$name"
    memory '50 GB'
    publishDir "${params.outdir}/qc/pileup", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_pileup
    file(bam_indices) from sorted_bam_indicies_for_pileup

    output:
    file("*.txt") into pileup_results

    script:
    """
    pileup.sh -Xmx20g \
              in=${bam_file} \
              out=${name}.coverage.stats.txt \
              hist=${name}.coverage.hist.txt
    """
 }

/*
 *STEP 5 - Create bedGraphs for visuzalization
 */

process bedgraphs {
    tag "$name"
    memory '80 GB'
    time '4h'
    publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.bedGraph"

    input:
    set val(name), file(bam_file) from sorted_bam_for_bedgraph
    set val(name), file(bam_indices) from sorted_bam_index_for_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph

    output:
    set val(name), file("${name}.bedGraph") into non_normalized_bedgraphs
    set val(name), file("${name}.rcc.bedGraph") into bedgraph_tdf, bedgraph_bigwig

    script:
    """
    genomeCoverageBed \
                     -bg \
                     -strand + \
                     -g ${params.chrom_sizes} \
                     -ibam ${bam_file} \
                     > ${name}.unsorted.bedGraph

    sortBed \
             -i ${name}.unsorted.bedGraph \
             > ${name}.bedGraph

    rm ${name}.unsorted.bedGraph

    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph \

    """
 }

/*
 *STEP 5+ - Convert normalized bedGraphs to bigWigs
 */

process normalized_bigwigs {
    tag "$name"
    memory '30 GB'
    publishDir "${params.outdir}/mapped/rcc_bigwig", mode: 'copy'
    
    when:
    params.savebw

    input:
    set val(name), file(bedgraph) from bedgraph_bigwig
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.rcc.bw") into normalized_bigwig

    script:
    """
    ${params.bedGraphToBigWig} ${bedgraph} ${chrom_sizes} ${name}.rcc.bw

    """
}

/*
 *STEP 5+ - IGV Tools : generate tdfs for optimal visualization in Integrative Genomics Viewer (IGV)
 */

process igvtools {
    tag "$name"
    memory '200 GB'
    time '1h'
    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/mapped/tdfs", mode: 'copy', pattern: "*.tdf"

    input:
    set val(name), file(normalized_bg) from bedgraph_tdf
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.tdf") into tiled_data_ch

    script:
    """
    /opt/igvtools/2.3.75/igvtools toTDF ${normalized_bg} ${name}.rcc.tdf ${chrom_sizes}
    """
 }



/*
 * STEP 6 - MultiQC
 */
process multiqc {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC

    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('qc/samtools_mapstats/*') from bam_flagstat.collect()
    file ('qc/rseqc/*') from rseqc_results.collect()
    file ('qc/preseq/*') from preseq_results.collect()
    file ('qc/hisat2_mapstats/*') from hisat2_mapstats.collect()
    file ('qc/picard/*') from picard_stats_multiqc.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}



/*
 * STEP 7 - Output Description HTML
 */
//
//process output_documentation {
//    tag "$prefix"
//    publishDir "${params.outdir}/pipeline_info/", mode: 'copy'
//
//    input:
//    file output_docs
//
//    output:
//    file "results_description.html"
//
//    script:
//    """
//    markdown_to_html.r $output_docs results_description.html
//    """
//}
//


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[ChIP-Flow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[ChIP-Flow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[ChIP-Flow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[ChIP-Flow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[ChIP-Flow] Pipeline Complete"

}
