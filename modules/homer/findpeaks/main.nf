process HOMER_FINDPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3" : null)
    
    publishDir "${params.outdir}/homer_find_peaks", mode: 'copy'

    input:
    tuple val(meta), path(tagDir)

    output:
    tuple val(meta), path("*peaks.txt"), emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    findPeaks \\
        $tagDir \\
        $args \\
        -o ${prefix}.peaks.txt
    """
}