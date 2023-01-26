/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def valid_params = [
//     ena_metadata_fields : ['run_accession', 'experiment_accession', 'library_layout', 'fastq_ftp', 'fastq_md5']
// ]

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// // Validate input parameters
// WorkflowSra.initialise(params, log, valid_params)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { FASTQC    } from '../modules/fastqc/main'
include { MULTIQC   } from '../modules/multiqc/main'

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
    RUN MAIN WORKFLOW
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
