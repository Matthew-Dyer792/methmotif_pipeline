process MULTIQC { 
    tag "pre-aligned_multiqc"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::multiqc=1.13' : null)
    container "${ workflow.containerEngine == 'singularity' ? 'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' : null}"

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    multiqc -f $args .
    """
}