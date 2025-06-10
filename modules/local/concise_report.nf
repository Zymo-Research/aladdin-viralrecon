process MULTIQC_PANGOLIN {
    label 'process_medium'

    conda "bioconda::multiqc=1.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.16--pyhdfd78af_0' :
        'biocontainers/multiqc:1.16--pyhdfd78af_0' }"

    input:
    path(multiqc_config)
    path(multiqc_logo)
    path('fastp/*')
    path('variants/*')
    path "multiqc_aladdin_viralrecon"

    output:
    path "*multiqc_report.html"     , emit: report
    path "*_data"                   , emit: data
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def title = params.project ? "--title \"Concise Aladdin Viralrecon Report for ${params.project}\"" : ''
    def filename = params.project ? "--filename " + params.project.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''

    """
    python -m venv venv
    source venv/bin/activate
    pip install -e multiqc_aladdin_viralrecon --no-cache-dir

    multiqc \\
        --force \\
        $title \\
        $filename \\
        $args \\
        $config \\
        .
    """
}