process FASTQC { 
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1' : null}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    if (meta.single_end) {
        """
        fastqc $args --threads $task.cpus ${reads}
        """
    } else {
        """
        fastqc $args --threads $task.cpus ${reads[0]} ${reads[1]}
        """
    }
}