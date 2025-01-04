process igv_reports {
    label 'high_cpu_high_plus_mem'
    container = params.containers.igv_reports
    tag "All Samples"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    path (bedgraphs)
    path (house_keeping_genes)
    path (genomeFile)
    path (genomeIndexFiles)
    tuple val(_), val(_),val(path_analysis),val(_), val(_)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 
    htmlFile = "igv_housekeeping_genes_report.html"

    output:
    path ("igv_housekeeping_genes_report.html")

    script:
    """
    create_report $house_keeping_genes --fasta $genomeFile --tracks $bedgraphs --output $htmlFile --flanking 100
    """
}