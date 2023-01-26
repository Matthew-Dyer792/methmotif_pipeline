process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::homer=4.11=pl526hc9558a2_3 bioconda::samtools=1.15" : null)

    publishDir "${params.outdir}/tag_dir", mode: 'copy'

    input:
    tuple val(meta), path(bed)
    path fasta

    output:
    tuple val(meta), path("${meta.id}"), emit: tagdir
    tuple val(meta), path("${meta.id}/tagInfo.txt"), emit: tag_info

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makeTagDirectory \\
        $prefix \\
        $args \\
        -genome $fasta \\
        $bed
    """
}