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

include { FASTQC as FASTQC_RAW_CHIP     } from '../modules/fastqc/main'
include { FASTQC as FASTQC_RAW_BS       } from '../modules/fastqc/main'
include { FASTQC as FASTQC_TRIM_CHIP    } from '../modules/fastqc/main'
include { FASTQC as FASTQC_TRIM_BS      } from '../modules/fastqc/main'
include { FASTP as FASTP_CHIP           } from '../modules/fastp/main'
include { STAR_ALIGN                    } from '../modules/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS       } from '../subworkflows/bam_sort_stats_samtools/main'
include { MULTIQC as MULTIQC_CHIP       } from '../modules/multiqc/main'
include { MULTIQC as MULTIQC_BS         } from '../modules/multiqc/main'

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// Create a channel for input read files
Channel
    .fromPath(params.chip_input, type: 'file', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .csv file matching: ${params.chip_input}\n" }
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
    .ifEmpty { exit 1, "Cannot find any .csv file matching: ${params.bs_input}\n" }
    .splitCsv(header: true)
    .map { row -> if (row.fastq_2) {
            tuple([id: "${row.sample}", single_end: false, trim: "${row.trim}".toBoolean()], [file("${row.fastq_1}"), file("${row.fastq_2}")]) 
        } else {
            tuple([id: "${row.sample}", single_end: true, trim: "${row.trim}".toBoolean()], file("${row.fastq_1}")) 
        } }
    .set { bs_reads }

// star index value channel
Channel
    .fromPath(params.star_index, type: 'dir', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any star index directory: ${params.star_index}\n" }
    .first()
    .set { star_index }

// samtools stats reference fasta value channel
Channel
    .fromPath(params.samtools_fasta, type: 'file', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any reference fasta file: ${params.samtools_fasta}\n" }
    .map { it -> tuple([id: it.simpleName], file(it)) }
    .first()
    .set { samtools_fasta }

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

//
// WORKFLOW: Run the pre-align steps
//
workflow ALIGN {
    //
    // STEP: Get fastqc for all chip reads
    //
    FASTQC_RAW_CHIP (chip_reads)
    
    // create a channel with reads to be trimmed
    chip_reads
        .filter { it.first().trim }
        .set{ chip_trim }

    //
    // STEP: trim chip reads where neccessary
    //
    FASTP_CHIP (chip_trim, params.fastp_adapter_fasta, params.fastp_save_trimmed_fail, params.fastp_save_merged)

    //
    // STEP: Recompute fastqc on trimmed chip reads
    //
    FASTQC_TRIM_CHIP (FASTP_CHIP.out.reads)

    // create a channel with trimmed and untrimmed reads to align
    chip_reads
        .filter { !it.first().trim }
        .mix( FASTP_CHIP.out.reads )
        .set{ prepared_chip_reads }
        // .view()

    //
    // STEP: Use star to align chip reads
    //
    STAR_ALIGN (prepared_chip_reads, star_index, params.star_gtf, params.star_ignore_sjdbgtf, params.star_seq_platform, params.star_seq_center)

    //
    // STEP: Use samtools to sort and index bam then produce stats
    //
    BAM_SORT_STATS_SAMTOOLS (STAR_ALIGN.out.bam_unsorted, samtools_fasta)

    // create a channel containing the multiqc files
    FASTQC_RAW_CHIP.out.zip
        .mix(FASTP_CHIP.out.json, FASTQC_TRIM_CHIP.out.zip, STAR_ALIGN.out.log_final)
        .map { it -> it[1] }
        .set{ mulitqc_chip_files }

    //
    // STEP: Summarize results with multiqc
    //
    MULTIQC_CHIP (mulitqc_chip_files, params.multiqc_config, params.extra_multiqc_config, params.multiqc_logo)

// ___________________________________________________________________________

    //
    // STEP: Get fastqc for all bsseq reads
    //
    // FASTQC_RAW_BS (bs_reads)

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
