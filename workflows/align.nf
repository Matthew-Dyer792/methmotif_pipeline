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

include { FASTQC as PRE_FASTQ_CHIP  } from '../modules/fastqc/main'
include { FASTQC as PRE_FASTQ_BS    } from '../modules/fastqc/main'
include { FASTQC as POST_FASTQ_CHIP } from '../modules/fastqc/main'
include { FASTQC as POST_FASTQ_BS   } from '../modules/fastqc/main'
include { FASTP as FASTP_CHIP       } from '../modules/fastp/main'
include { STAR_ALIGN                } from '../modules/star/align/main'
include { MULTIQC as MULTIQC_CHIP   } from '../modules/multiqc/main'
include { MULTIQC as MULTIQC_BS     } from '../modules/multiqc/main'

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// Create a channel for input read files
Channel
    .fromPath(params.chip_input, type: 'file', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .csv file matching: ${params.input}\n" }
    .splitCsv(header: true)
    .map { row -> if (row.fastq_2) {
            tuple([id: "${row.sample}", single_end: false, trim: "${row.trim}".toBoolean()], [file("${row.fastq_1}"), file("${row.fastq_2}")]) 
        } else {
            tuple([id: "${row.sample}", single_end: true, trim: "${row.trim}".toBoolean()], file("${row.fastq_1}")) 
        } }
    .set { chip_reads }

// Create a channel for input read files
Channel
    .fromPath(params.bs_input, type: 'file', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .csv file matching: ${params.input}\n" }
    .splitCsv(header: true)
    .map { row -> if (row.fastq_2) {
            tuple([id: "${row.sample}", single_end: false, trim: "${row.trim}".toBoolean()], [file("${row.fastq_1}"), file("${row.fastq_2}")]) 
        } else {
            tuple([id: "${row.sample}", single_end: true, trim: "${row.trim}".toBoolean()], file("${row.fastq_1}")) 
        } }
    .set { bs_reads }

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

//
// WORKFLOW: Run the pre-align steps
//
workflow ALIGN {
    // PRE_FASTQ_CHIP (chip_reads)
    
    chip_reads
        .filter { it.first().trim }
        // .set{ chip_trim }
        .view()

    // FASTP_CHIP (chip_trim)

    // POST_FASTQ_CHIP (FASTP_CHIP.out.reads)

    // filter out for multiqc... actually no

    chip_reads
        .filter { !it.first().trim }
        // .set{ chip_trim }
        .view()

    // STAR_ALIGN ()

    // MULTIQC_CHIP ()

    // PRE_FASTQ_BS (bs_reads)

    // FASTQC.out.zip
    //     .map{ it[1] }
    //     .collect()
    //     .set{ fastqc_zip }

    // fastqc_zip
    //     .set{ multiqc_files }

    // MULTIQC (multiqc_files)
}

/*
========================================================================================
    THE END
========================================================================================
*/
