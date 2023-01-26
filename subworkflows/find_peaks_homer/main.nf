//
// Finding peaks with Homer
//

params.tag_dir_options      = [:]
params.find_peaks__options  = [:]

include { HOMER_MAKETAGDIRECTORY    } from '../../modules/homer/maketagdirectory/main'  addParams( options: params.tag_dir_options )
include { HOMER_FINDPEAKS           } from '../../modules/homer/findpeaks/main'         addParams( options: params.find_peaks__options )

workflow FIND_PEAKS {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    reference_genome // channel: /path/to/reference/genome/.fa

    main:
    //
    // Build Tag Directories with Homer
    //
    HOMER_MAKETAGDIRECTORY ( bam, reference_genome )

    //
    // Find peaks in from the tag directories against the control (input)
    //
    HOMER_FINDPEAKS ( HOMER_MAKETAGDIRECTORY.out.tagdir )

    emit:
    tag_dir = HOMER_MAKETAGDIRECTORY.out.tagdir // channel: [ val(meta), bam   ]

    peaks   = HOMER_FINDPEAKS.out.txt           // channel: [ val(meta), [ bam ] ]
}