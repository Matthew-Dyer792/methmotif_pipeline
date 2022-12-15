#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    methmotif pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/fetchngs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WorkflowMain.initialise(workflow, params, log)

// // Check if --input file is empty
// ch_input = file(params.input, checkIfExists: true)
// if (ch_input.isEmpty()) {exit 1, "File provided with --input is empty: ${ch_input.getName()}!"}

// // Read in ids from --input file
// Channel
//     .from(file(params.input, checkIfExists: true))
//     .splitCsv(header:false, sep:'', strip:true)
//     .map { it[0] }
//     .unique()
//     .set { ch_ids }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// // Auto-detect input id type
// def input_type = ''
// if (WorkflowMain.isSraId(ch_input, log)) {
//     input_type = 'sra'
// } else if (WorkflowMain.isSynapseId(ch_input, log)) {
//     input_type = 'synapse'
// } else {
//     exit 1, 'Ids provided via --input not recognised please make sure they are either SRA / ENA / DDBJ or Synapse ids!'
// }

// if (params.input_type == input_type) {
//     if (params.input_type == 'sra') {
//         include { SRA } from './workflows/sra'
//     } else if (params.input_type == 'synapse') {
//         include { SYNAPSE } from './workflows/synapse'
//     }
// } else {
//     exit 1, "Ids auto-detected as ${input_type}. Please provide '--input_type ${input_type}' as a parameter to the pipeline!"
// }

if (params.workflow == 'pre-align') {
    include { PRE_ALIGN } from './workflows/pre-align'
} else if (params.workflow == 'trimming') {
    include { TRIMMING } from './workflows/trimming'
} else if (params.workflow == 'align') {
    include { TRIMMING } from './workflows/align'
}

//
// WORKFLOW: Run main methmotif pipeline depending on the step provided
//
workflow     METHMOTIF_PIPELINE {

    //
    // WORKFLOW: run the pre-alignment qc
    //
    if (params.workflow == 'pre-align') {
        PRE_ALIGN (  )

    //
    // WORKFLOW: trim the fastq files if necessary
    //
    } else if (params.workflow == 'trimming') {
        TRIMMING (  )

    //
    // WORKFLOW: trim the fastq files if necessary
    //
    } else if (params.workflow == 'align') {
        ALIGN (  )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    METHMOTIF_PIPELINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/