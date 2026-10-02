#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// input params
params.inputFile = "/home/achopra/BPA_Alt_Human/BPA/Fastq_dump/input/SRR_Acc_List.txt"
params.fd_outDir = "/home/achopra/BPA_Alt_Human/BPA/Fastq_dump/"
params.qc_outdir = "/home/achopra/BPA_Alt_Human/BPA/Fast_qc/"
params.qc_trimming = "/home/achopra/BPA_Alt_Human/BPA/Fast_p/"
params.multiqc = "/home/achopra/BPA_Alt_Human/BPA/Multi_qc/"
params.base_dir = "/home/achopra/BPA_Alt_Human/BPA/"
params.probes_dir = "/home/achopra/BPA_Alt_Human/BPA/probes/"
params.probes_output = "/home/achopra/BPA_Alt_Human/BPA/probes_output/"
params.sort_by_coordinate = "/home/achopra/BPA_Alt_Human/BPA/Picard_coordinate/"
params.mark_duplicate = "/home/achopra/BPA_Alt_Human/BPA/Mark_duplicates/"
params.picard_jar_path = "/home/achopra/miniconda3/envs/picard/share/picard-3.3.0-0/picard.jar"
params.feature_counts = "/home/achopra/BPA_Alt_Human/BPA/Feature_counts/"
params.annotation_path="/home/achopra/BPA_Alt_Human/BPA/Human_ref_files/temposeq_annotations_with_transcripts.gtf"

log.info """\
    =====================================================
       R N A S E Q - T E M P O S E Q   P I P E L I N E
    =====================================================
     READS : ${params.inputFile}
     """.stripIndent()

Channel
    .fromPath(params.inputFile)
    .splitText()
    .map { it.trim() }
    .set { FastqDump_ch }

process FastqDump {
    label 'Download_Raw_Reads'
    cpus 4
    clusterOptions "--nodes=1"
    tag "on ${sra_id}"
    publishDir "${params.fd_outDir}${sra_id}", mode: 'copy', pattern: "*.fastq"
    conda "/home/achopra/miniconda3/envs/fastqdump/fastqdump.yaml"                

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}.fastq")

    script:
    """
    echo "Downloading: ${sra_id}" 
    prefetch ${sra_id}
    fasterq-dump ${sra_id} -e 4 --skip-technical
    """
}

process QualityCheck {
    label "Quality_Checks" 
    cpus 4
    clusterOptions "--nodes=1" 
    tag "on ${sra_id}"
    publishDir "${params.qc_outdir}${sra_id}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/fastqc/fastqc.yaml"

    input:
    tuple val(sra_id), path(fastq_file)

    output:
    tuple val(sra_id), path("${sra_id}_fastqc.html"), path("${sra_id}_fastqc.zip")

    script:
    """
    echo "Processing fastqc file: ${sra_id}" 
    fastqc ${fastq_file} --outdir . --threads 4
    """
}

process FastpTrimming {
    label "Trimming_Adapter" 
    cpus 8
    clusterOptions "--nodes=1" 
    tag "on ${sra_id}"
    publishDir "${params.qc_trimming}${sra_id}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/fastp/fastp.yaml"

    input:
    tuple val(sra_id), path(fastq_file)

    output:
    tuple  val(sra_id), path("${sra_id}_fastp_trimmed.fastq"), path("${sra_id}_fastp.html"), path("${sra_id}_fastp.json")

    script:
    """
    fastp \
        -i ${fastq_file} \
        -o ${sra_id}_fastp_trimmed.fastq \
        -h ${sra_id}_fastp.html \
        -j ${sra_id}_fastp.json \
        --thread 8 \
        --cut_tail \
        --cut_tail_window_size 2 \
        --cut_tail_mean_quality 30 \
        --length_required 30
    """
}

process MultiQCFastqcData {
    label "MultiQC_On_QCfiles" 
    cpus 4
    clusterOptions "--nodes=1" 
    publishDir "${params.multiqc}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/multiqc/multiqc.yaml"

    input:
    path(zip_files)

    output:
    path("TESTmultiqc_FastQC_report.html")

    script:
    """
    multiqc ${zip_files} -n TESTmultiqc_FastQC_report.html -o . 
    """
}
process MultiQCFastPData {
    label "MultiQC_On_Pfiles"
    cpus 4
    clusterOptions "--nodes=1"
    publishDir "${params.multiqc}", mode: 'copy'
    conda "/home/achopra/miniconda3/envs/multiqc/multiqc.yaml"

    input:
    path(fastp_files)

    output:
    path("TESTfastp_multiqc_report.html")

    script:
    """
    multiqc ${fastp_files.join(' ')} -n TESTfastp_multiqc_report.html -o .
    """
}

process RNABowtie {
    label "Bowtie_probe_Aligning"
    cpus 28
    clusterOptions "--nodes=1"
    publishDir "${params.probes_output}${sra_id}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/bowtie/bowtie.yaml"
    memory '60 GB'

    input:
    tuple val(sra_id), path(trimmed_fastq), path(html_file), path(json_file)

    output:
    tuple val(sra_id), path("${sra_id}.sam"), path("${sra_id}_bowtie_Aligned.out.bam"), path("${sra_id}bowtie.log")

    script:
    """
    mkdir -p ${params.probes_output}${sra_id}
    bowtie2 --no-unal \
            -x ${params.probes_dir}probe_index \
            -U ${trimmed_fastq} \
            -S ${sra_id}.sam \
            -p ${task.cpus} > ${sra_id}bowtie.log 2>&1
    samtools view -bS ${sra_id}.sam | samtools sort -o ${sra_id}_bowtie_Aligned.out.bam
    # rm -f ${sra_id}.sam
    """
}

process MultiQCAlignedData {
    label "MultiQC_On_Alignedfiles" 
    cpus 4
    clusterOptions "--nodes=1" 
    publishDir "${params.multiqc}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/multiqc/multiqc.yaml"

    input:
    path(log_files)

    output:
    path("TESTmultiqc_Aligned_report.html")

    script:
    """
    multiqc ${log_files} -n TESTmultiqc_Aligned_report.html -o . 
    """
}

process SortByCoordinate {
    label "Sort_By_Coordinate"
    cpus 8
    clusterOptions "--nodes=1"
    tag "on ${sra_id}"
    publishDir "${params.sort_by_coordinate}${sra_id}", mode:'copy'
    conda "/home/achopra/miniconda3/envs/picard/picard.yaml"

    input:
    tuple val(sra_id), path(sam), path(aligned_bam), path(log_file)

    output:
    tuple val(sra_id), path("${sra_id}SC.bam")

    script:
    """
    java -jar ${params.picard_jar_path} SortSam \
        INPUT=${aligned_bam} \
        OUTPUT=${sra_id}SC.bam \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT 
    """
}

process FeatureCounts {
    label "Feature_Counts"
    cpus 8
    clusterOptions "--nodes=1"
    publishDir "${params.feature_counts}", mode:'move'
    conda "/home/achopra/miniconda3/envs/featurecounts/featurecounts.yaml"

    input:
    path(bam_files)

    output:
    path("TESTcombined_feature_counts.txt")

    script:
    """
    featureCounts -T 8 \
              -t exon \
              -g gene_id \
              -a "${params.annotation_path}" \
              -o "TESTcombined_feature_counts.txt" ${bam_files} \
              2> featurecounts.log
    """
}

workflow {
    FastqDump_Files = FastqDump(FastqDump_ch)
    
    Fastq_QC_Files = QualityCheck(FastqDump_Files)
    
    Trimmed_files = FastpTrimming(FastqDump_Files)
    
    QC_Zip_Files = Fastq_QC_Files.map { sra_id, html_file, zip_file -> zip_file }
    MultiQC_file_qc_data = MultiQCFastqcData(QC_Zip_Files.collect())
    
    P_Zip_Files = Trimmed_files.map { sra_id, trimmed_fastq, fastp_html, fastp_json -> [fastp_html, fastp_json] }.flatten()
MultiQC_file_P_data = MultiQCFastPData(P_Zip_Files.collect())


    
    Aligned_files = RNABowtie(Trimmed_files)
    
    Aligned_Log_Files = Aligned_files.map { sra_id, sam_file, aligned_bam, log_file -> log_file }
    MultiQC_file_aligned_data = MultiQCAlignedData(Aligned_Log_Files.collect())
    
    Sort_by_coordinate_files = SortByCoordinate(Aligned_files)
    
    bam_files = Sort_by_coordinate_files.map { sra_id, sorted_bam -> sorted_bam }
    FeatureCounts(bam_files.collect())
    
    workflow.onComplete {
        println("Pipeline completed successfully.")
    }
}