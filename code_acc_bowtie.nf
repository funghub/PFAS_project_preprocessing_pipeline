#!/usr/bin/env nextflow

// change file name to code.nf when finished
// REMEMBER TO CD into PFAS_Data_NF

// Command used, but now run through github, replace file with link
// nextflow run funghub/PFAS_project --input PRJNA604830 -profile spartan_hpc -latest -resume
// nextflow run code.nf --input PRJNA604830 -profile spartan_hpc -resume // this is for code stored on HPC

// nextflow.enable.dsl=2

process header {
    script:
    """
    datenow=\$(date)
    echo \$datenow
    # srun hostname # print the name of the node

    start=\$(date +%s)
    echo "start time: \$start" # Print the timestamp in seconds
    echo "hostname: \$HOSTNAME"
    echo ""

    """
}

process footer {
    script:
    """
    end=\$(date +%s)
    echo "end time: \$end" # Print the timestamp in seconds
    runtime=\$((end - start))
    echo "run time: \$runtime" # print the run time
    """
}

process retrieve_accessions_numbers {
    conda "bioconda::entrez-direct"

    input: 
    val SRA_accession_number

    output:
    path "SRR_Acc_List.txt", emit: accession_numbers_file

    script:
    """
    esearch -db sra -query ${SRA_accession_number} | efetch -format runinfo | cut -d',' -f 1 | grep SRR > SRR_Acc_List.txt
    """
    // you can do grep -c to double check on the web if the number of items match
}

process retrieve_sra {
    conda "bioconda::sra-tools=3.4.1 conda-forge::ossuuid"

    input:
    val accession_number

    output:
    path "${accession_number}/${accession_number}.sra", emit: sra_files


    script:
    """
    # download all SRA accessions
    # prefetch --option-file ${accession_number}
    prefetch ${accession_number} --output-directory .
    """
}

process retrieve_fastq {
    conda "bioconda::sra-tools=3.4.1 conda-forge::ossuuid"

    input:
    path sra_file // path not val because receiving .sra file from channel

    output:
    path "*.fastq", emit: pretrim_fastq

    script:
    """    
    # read each line from txt file and into dump 1 at a time (slow need to separate file)
    # convert SRA to FASTQ file one at a time
    fasterq-dump ${sra_file} --threads ${task.cpus}
    """
}

process fastqc_pretrim {
    conda "bioconda::fastqc=0.12.1" // conda-forge::perl=5.32.0"    
    
    input:
    path pretrim_fastq

    output:
    // path "*_fastqc.{zip,html}", emit: pretrim_qc_files
    path "*.{zip,html}", emit: pretrim_qc_files

    script:
    """
    fastqc -t ${task.cpus} ${pretrim_fastq}
    """
}

process MULTIQC_pretrim {
    conda "bioconda::multiqc"
    
    input:
    path pretrim_qc_files

    output:
    path "multiqc_report.html", emit: report_pretrim

    script:
    """
    multiqc ${pretrim_qc_files}
    """
}

process FASTP {
    conda "bioconda::fastp"
    // publishDir "PFAS_Data_NF/${task.index}/fastp_results", mode: 'copy'

    input:
    path fastq

    output:
    path "${fastq.baseName}_trimmed.fastq", emit: trimmed
    path "*.{json,html}"                  , emit: reports

    script:
    """
    fastp -i ${fastq} -o ${fastq.baseName}_trimmed.fastq \
          -j ${fastq.baseName}.fastp.json \
          -h ${fastq.baseName}.fastp.html
    """
}

process FASTQC {
    conda "bioconda::fastqc=0.12.1" // using latest version might rid error of java.lang openjdk
    // publishDir "PFAS_Data_NF/${task.index}/fastqc_results", mode: 'copy'

    input:
    path trimmed_fastq

    output:
    path "*_fastqc.{zip,html}", emit: qc_files

    script:
    """
    fastqc -t ${task.cpus} ${trimmed_fastq}
    """
}

process MULTIQC {
    conda "bioconda::multiqc"
    // publishDir "PFAS_Data_NF/${task.index}/multi_qc_results", mode: 'copy'
    
    input:
    // path 'stats/*'
    path reports

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc ${reports}
    """
}

process STAR_index {
    conda "bioconda::star"
    
    output:
    path "STAR_hg38_index", emit: star_index
    path "Genome_Annotation/hg38.ncbiRefSeq.gtf", emit: gtf_file

    script:
    """
    # get genome assembly and get only chr files
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
    mkdir Genome_Assembly
    tar -xvf hg38.chromFa.tar.gz -C Genome_Assembly
    ls Genome_Assembly/chroms | grep -v -e "random" -e "alt" -e "chrUn" | xargs -I{} cat Genome_Assembly/chroms/{} > Genome_Assembly/chroms_all.fa

    # get gene annotation gtf file
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
    mkdir Genome_Annotation
    gunzip -c hg38.ncbiRefSeq.gtf.gz > Genome_Annotation/hg38.ncbiRefSeq.gtf

    mkdir STAR_hg38_index
    
    # create index
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir STAR_hg38_index \
        --genomeFastaFiles Genome_Assembly/chroms_all.fa \
        --sjdbGTFfile Genome_Annotation/hg38.ncbiRefSeq.gtf \
        --sjdbOverhang 99
    """
}

process STAR_align {
    conda "bioconda::star"
    
    input:
    path star_index
    path trimmed_fastq

    output:
    path "*.bam", emit: star_alignment
    path "*.{out,tab}", emit: star_logs

    script:
    """
    # complete alignment
    STAR --genomeDir ${star_index} \
        --readFilesIn ${trimmed_fastq} \
        --outFileNamePrefix ${trimmed_fastq.baseName}. \
        --runThreadN ${task.cpus} \
        --outSAMtype BAM SortedByCoordinate
    """
}

process samtools_index { // for getting bai file from bam
    conda "bioconda::samtools"
    
    input:
    path star_alignment

    output:
    // path "${star_alignment.baseName}.bai", emit: bai_files
    path "*.bai", emit: bai_files // for a BAM file aln.bam, either aln.bam.bai will be created

    script:
    """
    samtools index ${star_alignment}
    """
}

process samtools_flagstat {
    conda "bioconda::samtools"
    
    input:
    path star_alignment

    output:
    path "*.txt", emit: flagstat

    script:
    def prefix = star_alignment.name.replace("_trimmed.Aligned.sortedByCoord.out.bam","")

    """
    samtools flagstat ${star_alignment} > ${prefix}_flagstat.txt
    """
}


process picard_add_read_groups {
    conda "bioconda::picard=3.3.0"
    
    input:
    path star_alignment

    output:
    path "*.bam", emit: add_RG_bam

    script:
    def prefix = star_alignment.name.replace("_trimmed.Aligned.sortedByCoord.out.bam","")

    """
    # make sure to add back the read groups to the header of the BAM file
    picard AddOrReplaceReadGroups I=${star_alignment} O=${prefix}_RG.bam RGID=${prefix} RGLB=lib_${prefix} RGPL=ILLUMINA RGPU=${prefix} RGSM=${prefix}
    """
}

process picard_mark_duplicates {
    conda "bioconda::picard=3.3.0"
    
    input:
    path add_RG_bam

    output:
    path "*.bam", emit: marked_dups_bam
    path "*metrics.txt", emit: marked_dups_metrics

    script:
    def prefix = add_RG_bam.name.replace("_RG.bam","")

    """
    picard MarkDuplicates I=${add_RG_bam} O=${prefix}_duplicates.bam M=${prefix}_dup_metrics.txt
    """
}

process MULTIQC_markdups_flagstat {
    conda "conda-forge::polars-lts-cpu bioconda::multiqc=1.33"
    
    input:
    path marked_dups_flagstat_metrics

    output:
    path "multiqc_markdups_flagstat.html", emit: report_markdups_flagstat

    script:
    """
    multiqc ${marked_dups_flagstat_metrics} -n multiqc_markdups_flagstat
    """
}

process feature_counts_raw {
    conda "bioconda::subread"
    
    input:
    path bam_files
    path gtf_file
    val prefix      // because of Module aliases, add prefix
    // val prefix is for _raw or _markdups (added as an input to when calling the process in workflow)

    output:
    path "${prefix}_counts.txt", emit: counts
    path "${prefix}_counts.txt.summary", emit: summary

    script:
    """
    echo hello # remove
    featureCounts -T ${task.cpus} -a ${gtf_file} -t exon -g gene_id -o ${prefix}_counts.txt ${bam_files}
    """
}

process feature_counts_markdups {
    conda "bioconda::subread"
    
    input:
    path bam_files
    path gtf_file
    val prefix      // because of Module aliases, add prefix, forget it: Module aliases must be from a separate .nf
    // val prefix is for _raw or _markdups (added as an input to when calling the process in workflow)

    output:
    path "${prefix}_counts.txt", emit: counts
    path "${prefix}_counts.txt.summary", emit: summary

    script:
    """
    echo hello # remove
    featureCounts -T ${task.cpus} -a ${gtf_file} -t exon -g gene_id -o ${prefix}_counts.txt ${bam_files}
    """
}


/*
 * Pipeline parameters
 */
// params {
//     // input: Path = '/scratch/home/lfung/PFAS_Data_NF/test_multiqc/*.fastq'

//     // input: Path = 'test_multiqc/*.fastq' // assuming cd into test_runs
//     input = 'test_multiqc/*.fastq' // assuming cd into test_runs (don't need path b/c already channel.fromPath)

// }

params.output_dir = 'results'     // Default output directory
// params.input = "${params.input_dir}/*.fastq"


params.help = false  // Set the default to false
include { paramsHelp } from 'plugin/nf-schema'

// Module aliases to reuse processes but have different outputs!: didn't work
// feature_counts must be located in separate .nf file
// https://training.nextflow.io/2.0/basic_training/modules/#module-aliases
// https://stackoverflow.com/questions/76730547/how-to-reuse-the-same-process-twice-in-within-the-same-module-in-nextflow-dsl2
// include { feature_counts as feature_counts_raw } from './modules/feature_counts'
// include { feature_counts as feature_counts_markdups } from './modules/feature_counts'


workflow {

    main:

    if (params.help) {
    log.info """
    FASTQ files -> FASTP trimmed files -> FASTQC files-> MULTIQC report Pipeline
    ----------------------
    Usage:
    nextflow run funghub/PFAS_project --input [SRA accession number] --output_dir [folder for results] -profile spartan_hpc -latest -resume
    
    Usage (default):
    nextflow run funghub/PFAS_project --input PRJNA604830 -profile spartan_hpc -latest -resume

    Options:
      --input    Input SRA accession number (no default value)
      --output_dir   Directory to save results (default: results)
    """
    exit 0
    }

    def sra_accession_number = params.input

    // header()

    retrieve_accessions_numbers(sra_accession_number)

    // split the txt file into one accession number per channel
    // this replaces xargs in the bash script in retrieve_fastq which was inefficient
    retrieve_accessions_numbers.out.accession_numbers_file
        .splitText()
        .map { it.trim() } // for each item in txt file, trim off white spaces char
        .filter {it} // get lines where it (each item) is true
        .set { accessions_ch } // set it as a new variable for channel

    retrieve_sra(accessions_ch)
    retrieve_fastq(retrieve_sra.out.sra_files)
    fastqc_pretrim(retrieve_fastq.out.pretrim_fastq)
    MULTIQC_pretrim(fastqc_pretrim.out.pretrim_qc_files.collect())
    
    // comment out below if need to verify pretrim multiqc report first 
    FASTP(retrieve_fastq.out.pretrim_fastq)
    FASTQC(FASTP.out.trimmed) // .trimmed specifically refers to the emit name given
    MULTIQC(FASTQC.out.qc_files.collect()) // must use .colect() with () to work
    STAR_index()
    STAR_align(STAR_index.out.star_index, FASTP.out.trimmed)

    samtools_index(STAR_align.out.star_alignment)
    samtools_flagstat(STAR_align.out.star_alignment)

    picard_add_read_groups(STAR_align.out.star_alignment)
    picard_mark_duplicates(picard_add_read_groups.out.add_RG_bam)

    // add in picard metrics file and mix channel with the outputs for samtools flagstat metrics
    MULTIQC_markdups_flagstat(picard_mark_duplicates.out.marked_dups_metrics.collect().mix(samtools_flagstat.out.flagstat.collect()).collect())
    
    // feature counts for without marked duplications!!!
    feature_counts_raw(STAR_align.out.star_alignment.collect(), STAR_index.out.gtf_file, "raw")
    // feature counts for with marked duplications!!!
    feature_counts_markdups(picard_mark_duplicates.out.marked_dups_bam.collect(), STAR_index.out.gtf_file, "markdups")
    


    // footer()


    publish:

    retrieve_accessions_numbers = retrieve_accessions_numbers.out.accession_numbers_file

    sra_files = retrieve_sra.out.sra_files
    fastq_pretrim = retrieve_fastq.out.pretrim_fastq
    fastqc_results_pretrim = fastqc_pretrim.out.pretrim_qc_files
    multiqc_results_pretrim = MULTIQC_pretrim.out.report_pretrim

    // comment out below if need to verify pretrim multiqc report first 
    fastp_trimmed = FASTP.out.trimmed
    fastp_reports = FASTP.out.reports

    fastqc_results = FASTQC.out.qc_files
    
    multiqc_results = MULTIQC.out.report

    star_index = STAR_index.out.star_index
    star_alignment = STAR_align.out.star_alignment
    star_logs = STAR_align.out.star_logs

    bai_files = samtools_index.out.bai_files

    flagstat = samtools_flagstat.out.flagstat

    add_read_groups_bam = picard_add_read_groups.out.add_RG_bam

    marked_dups_bam = picard_mark_duplicates.out.marked_dups_bam
    marked_dups_metrics = picard_mark_duplicates.out.marked_dups_metrics

    multiqc_markdups_flagstat = MULTIQC_markdups_flagstat.out.report_markdups_flagstat

    // feature counts for without marked duplications!!!
    featurecounts_raw = feature_counts_raw.out.counts
    featurecounts_summary_raw = feature_counts_raw.out.summary

    featurecounts_markdups = feature_counts_markdups.out.counts
    featurecounts_summary_markdups = feature_counts_markdups.out.summary

}

output {

    retrieve_accessions_numbers {
        path "${params.output_dir}/pretrim/SRR_Acc_List.txt"
    }

    sra_files {
        path "${params.output_dir}/pretrim/sra_files"
    }

    fastq_pretrim {
        path "${params.output_dir}/pretrim/fastq_pretrim" // tells where to save outputs to
    }

    fastqc_results_pretrim {
        path "${params.output_dir}/pretrim/fastqc_results_pretrim"
    }

    multiqc_results_pretrim {
        path "${params.output_dir}/pretrim/multi_qc_results_pretrim"
        mode 'copy'
    }

    // comment out below if need to verify pretrim multiqc report first 
    // KEEP
    fastp_trimmed {
        path "${params.output_dir}/fastp_results/trimmed" // tells where to save outputs to
        // without mode copy below, we are softlinking
        mode 'copy'
    }

    fastp_reports {
        path "${params.output_dir}/fastp_results/reports" // tells where to save outputs to
        // without mode copy below, we are softlinking
        // mode 'copy'
    }

    // KEEP
    fastqc_results {
        path "${params.output_dir}/fastqc_results"
        // mode 'copy'
    }

    multiqc_results {
        path "${params.output_dir}/multi_qc_results"
        mode 'copy'
    }

    star_index {
        path "${params.output_dir}/STAR_hg38_index"
        // mode 'copy'
    }

    star_alignment {
        path "${params.output_dir}/STAR_alignment"
        mode 'copy'
    }

    star_logs {
        path "${params.output_dir}/STAR_logs"
        // mode 'copy'
    }

    // all the BAI files per BAM file
    bai_files {
        path "${params.output_dir}/STAR_alignment"
        mode 'copy'
    }

    flagstat {
        path "${params.output_dir}/STAR_flagstat"
        // mode 'copy'
    }

    add_read_groups_bam {
        path "${params.output_dir}/picard/add_RG_bam"
        mode 'copy'
    }

    marked_dups_bam {
        path "${params.output_dir}/picard/marked_dups_bam"
        mode 'copy'
    }

    marked_dups_metrics {
        path "${params.output_dir}/picard/marked_dups_metrics"
        mode 'copy'
    }

    multiqc_markdups_flagstat {
        path "${params.output_dir}/multi_qc_results"
        mode 'copy'
    }

    featurecounts_raw {
        path "${params.output_dir}/featurecounts"
        mode 'copy'
    }

    featurecounts_summary_raw {
        path "${params.output_dir}/featurecounts"
        mode 'copy'
    }

    featurecounts_markdups {
        path "${params.output_dir}/featurecounts"
        mode 'copy'
    }

    featurecounts_summary_markdups {
        path "${params.output_dir}/featurecounts"
        mode 'copy'
    }

}
