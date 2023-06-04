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
include { FASTP as FASTP_BS             } from '../modules/fastp/main'
include { STAR_ALIGN                    } from '../modules/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS       } from '../subworkflows/bam_sort_stats_samtools/main'
include { BISMARK_ALIGN                 } from '../modules/bismark/align/main'
include { BISMARK_DEDUPLICATE           } from '../modules/bismark/deduplicate/main'
include { BISMARK_METHYLATIONEXTRACTOR  } from '../modules/bismark/methylationextractor/main'
include { MULTIQC as MULTIQC_CHIP       } from '../modules/multiqc/main'
include { MULTIQC as MULTIQC_BS         } from '../modules/multiqc/main'

/*
========================================================================================
    BUILD WORKFLOW CHANNELS
========================================================================================
*/

// Create a channel for input chip read files
Channel
    .fromPath(params.chip_input, type: 'file', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any .csv file matching: ${params.chip_input}\n" }
    .splitCsv(header: true)
    .map { row -> if (row.fastq_2) {
            tuple([id: "${row.sample}", single_end: false, trim: "${row.trim}".toBoolean(), control: "${row.control}".toBoolean()], [file("${row.fastq_1}"), file("${row.fastq_2}")]) 
        } else {
            tuple([id: "${row.sample}", single_end: true, trim: "${row.trim}".toBoolean(), control: "${row.control}".toBoolean()], file("${row.fastq_1}")) 
        } }
    .set { chip_reads }

// Create a channel for input bsseq read files
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

// bsimark index value channel
Channel
    .fromPath(params.bismark_index, type: 'dir', checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any bismark index directory: ${params.bismark_index}\n" }
    .first()
    .set { bismark_index }

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
    BAM_SORT_STATS_SAMTOOLS (STAR_ALIGN.out.bam, samtools_fasta)

    // create a channel containing the multiqc files
    FASTQC_RAW_CHIP.out.zip
        .mix(FASTP_CHIP.out.json, FASTQC_TRIM_CHIP.out.zip, STAR_ALIGN.out.log_final)
        .map { it -> it[1] }
        .collect()
        .set{ mulitqc_chip_files }

    //
    // STEP: Summarize results with multiqc
    //
    MULTIQC_CHIP (mulitqc_chip_files, params.multiqc_config, params.extra_multiqc_config, params.multiqc_logo)

    // create a channel to split ids by an "_"
    STAR_ALIGN.out.bam
        .map{ it -> tuple([id: it[0].id.split("_")[0], single_end: it[0].single_end, trim: it[0].trim, control: it[0].control], it[1]) }
        .set{ adjusted_ids }

    // create a channel with alignments that are experimental chip-seq
    adjusted_ids
        .filter { !it.first().control }
        .groupTuple()
        .set{ tf_chip }

    // create a channel with alignments that are control chip-seq
    adjusted_ids
        .filter { it.first().control }
        .groupTuple()
        .first()
        .set{ control_chip }

    //
    // STEP: Encode-dcc pipeline for 
    //
    MULTIQC_CHIP (tf_chip, control_chip)

// ___________________________________________________________________________

    //
    // STEP: Get fastqc for all bsseq reads
    //
    FASTQC_RAW_BS (bs_reads)

    // create a channel with reads to be trimmed
    bs_reads
        .filter { it.first().trim }
        .set{ bs_trim }

    //
    // STEP: trim bsseq reads where neccessary
    //
    FASTP_BS (bs_trim, params.fastp_adapter_fasta, params.fastp_save_trimmed_fail, params.fastp_save_merged)

    //
    // STEP: Recompute fastqc on trimmed bsseq reads
    //
    FASTQC_TRIM_BS (FASTP_BS.out.reads)

    // create a channel with trimmed and untrimmed reads to align
    bs_reads
        .filter { !it.first().trim }
        .mix( FASTP_BS.out.reads )
        .set{ prepared_bs_reads }

    //
    // STEP: Use bismark (bowtie2) to align bsseq reads
    //
    BISMARK_ALIGN (prepared_bs_reads, bismark_index)

    //
    // STEP: Use bismark
    //
    BISMARK_DEDUPLICATE (BISMARK_ALIGN.out.bam)

    //
    // STEP: Use bismark
    //
    BISMARK_METHYLATIONEXTRACTOR (BISMARK_DEDUPLICATE.out.bam, bismark_index)

    // create a channel containing the multiqc files
    FASTQC_RAW_BS.out.zip
        .mix(FASTP_BS.out.json, FASTQC_TRIM_BS.out.zip, BISMARK_ALIGN.out.report, BISMARK_DEDUPLICATE.out.report, BISMARK_METHYLATIONEXTRACTOR.out.report, BISMARK_METHYLATIONEXTRACTOR.out.mbias)
        .map { it -> it[1] }
        .collect()
        .set{ mulitqc_bsseq_files }

    //
    // STEP: Summarize results with multiqc
    //
    MULTIQC_BS (mulitqc_bsseq_files, params.multiqc_config, params.extra_multiqc_config, params.multiqc_logo)
}

/*
========================================================================================
    THE END
========================================================================================
*/
