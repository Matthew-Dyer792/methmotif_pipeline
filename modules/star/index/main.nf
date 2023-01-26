process STAR_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::star=2.7.10b" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/star:2.7.10b--h9ee0642_0' : null }"

    input:
    path fasta

    output:
    path "star", emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_list = args.tokenize()
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        --runThreadN $task.cpus \\
        $memory \\
        $args
    """
}