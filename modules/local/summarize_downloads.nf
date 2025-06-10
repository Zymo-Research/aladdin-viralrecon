process SUMMARIZE_DOWNLOADS {
    cache false
    publishDir "${params.outdir}/download_data", mode: 'copy'

    input:
    path locations
    path design

    output:
    path 'files_to_download.json'

    script:
    """
    summarize_downloads.py $locations -d $design
    """
}