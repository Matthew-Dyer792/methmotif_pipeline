process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::star=2.7.10b" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/star:2.7.10b--h9ee0642_0' : null }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_sam_type = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    def mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --runThreadN $task.cpus \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $args
    $mv_unsorted_bam
    """
}