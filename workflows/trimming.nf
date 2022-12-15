#!/usr/bin/env nextflow
/*
========================================================================================
    methmotif trimming
========================================================================================
    remove the stated sequences from the provided fastq files. Then re-run quality 
    control to confirm the changes were made correctly. Multiqc report provided.
----------------------------------------------------------------------------------------
*/


//needs rework




nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE
========================================================================================
*/

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// Create a channel for input read files
Channel
    .fromPath(params.input, type: 'file', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .csv file matching: ${params.input}\n" }
    .splitCsv(header: true)
    .map { row -> if (row.fastq_2) {
            tuple([id: "${row.sample}", single_end: false, phred: "${row.phred}", trim: row.trim == 'yes' ? true : false], [file("${row.fastq_1}"), file("${row.fastq_2}")]) 
        } else {
            tuple([id: "${row.sample}", single_end: true, phred: "${row.phred}", trim: row.trim == 'yes' ? true : false], file("${row.fastq_1}")) 
        } }
    .set { in_file }

// institute a filtering step to remove the files that do not need to be trimmed
in_file
    .filter { it[0].trim == true }
    .set { reads }
/*
========================================================================================
    PIPELINE STEPS
========================================================================================
*/

process TRIMMOMATIC { 
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' : null}"

    // publishDir "${params.outdir}/trimmed_fastq", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: trimmed_fastq
    path('*_unpaired.fastq.gz'), optional:true, emit: unpaired_fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    if (meta.single_end) {
        """
        trimmomatic SE \\
            -threads $task.cpus \\
            ${meta.phred} \\
            ${reads} \\
            ${prefix}_trimmed.fastq.gz \\
            $args
        """
    } else {
        """
        trimmomatic PE \\
            -threads $task.cpus \\
            ${meta.phred} \\
            ${reads[0]} ${reads[1]} \\
            ${prefix}_1_trimmed.fastq.gz ${prefix}_1_trimmed_unpaired.fastq.gz \\
            ${prefix}_2_trimmed.fastq.gz ${prefix}_2_trimmed_unpaired.fastq.gz\\
            $args
        """
    }
}

//possibly add cutadapt or trim-galore to work trim WGBS files

//work on moving fastqc and multiqc into modules to run shared between pre-align and trimming

process FASTQC { 
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1' : null}"

    // publishDir "${params.outdir}/pre-align_fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    if (meta.single_end) {
        """
        fastqc $args --threads $task.cpus ${reads}
        """
    } else {
        """
        fastqc $args --threads $task.cpus ${reads[0]} ${reads[1]}
        """
    }
}

process MULTIQC { 
    tag "pre-aligned_multiqc"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.13' : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' : null}"

    // publishDir "${params.outdir}/pre-align_multiqc", mode: 'copy'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    multiqc -f $args .
    """
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Run the pre-align steps
//
workflow TRIMMING {
    TRIMMOMATIC (reads)

    FASTQC (TRIMMOMATIC.out.trimmed_fastq)

    //integrate the trimmed data back in with the non-trimmed data for a global view

    FASTQC.out.zip
        .map{ it[1] }
        .collect()
        .set{ fastqc_zip }

    fastqc_zip
        .set{ multiqc_files }

    MULTIQC (multiqc_files)
}

/*
========================================================================================
    THE END
========================================================================================
*/
