#!/usr/bin/env nextflow
/*
========================================================================================
    methmotif pre-align
========================================================================================
    take in or downnload fastq files and provide out a multiqc report on the data
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
            tuple([id: "${row.sample}", single_end: false], [file("${row.fastq_1}"), file("${row.fastq_2}")]) 
        } else {
            tuple([id: "${row.sample}", single_end: true], file("${row.fastq_1}")) 
        } }
    .set { reads }

/*
========================================================================================
    PIPELINE STEPS
========================================================================================
*/

// possible fetch-ngs step

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
workflow PRE_ALIGN {
    FASTQC (reads)

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
