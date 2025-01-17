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

process igv_sample_reports {
    label 'high_cpu_high_plus_mem'
    container = params.containers.igv_reports
    tag "Sample - $sampleId"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    tuple val(sampleId),val(path_analysis),path(bedgraph),val (_)
    each path (house_keeping_genes)
    each path (genomeFile)
    each path (genomeIndexFiles)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 
    htmlFile = sampleId + "_igv_housekeeping_genes_report.html"

    output:
    path ("*.html")

    script:
    """
    create_report $house_keeping_genes --fasta $genomeFile --tracks $bedgraph --output $htmlFile 
    """
}

process igv_consolidate_report {
    label 'low_cpu_low_plus_mem'
    container = params.containers.ubuntu
    tag "All Samples"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    tuple val(_),val(_),val(path_analysis),val(_),val (_)
    path (samples_report)
    path (house_keeping_header)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 
    htmlFile = "igv_housekeeping_genes_mqc.html"

    output:
    path ("*.html")

    script:
    """
    # Define the output HTML file
    cat $house_keeping_header > $htmlFile

    # Append HTML links for matching files to the consolidated HTML file
    for file in *_igv_housekeeping_genes_report.html; do
        # Extract the string before '_igv_housekeeping_genes_report.html'
        link_text=\$(basename "\$file" "_igv_housekeeping_genes_report.html")
        echo "<a href='\${file}' target='_blank' class='btn btn-primary'>\${link_text}</a>" >> $htmlFile
    done
    """
}

process igv_session {
    label 'med_cpu_med_mem'
    container = params.containers.python
    tag "All Samples"

    publishDir "$path_sample_igv", mode : 'copy'

    input:
    tuple val(_),val(_),val(path_analysis),val(_),val (_)
    path (tdf)
    path (chIGVFilestoSessions)
    tuple val(genome), val(_), val(_), val(_), val(_)
    path (house_keeping_genes)

    exec:
    path_sample_igv = path_analysis + "/reports/igv_session/"
    fileOut = "IGV_Session.xml"

    output:
    tuple path ("*.xml"),path (bedgraph)

    script:
    """
    python $chIGVFilestoSessions $fileOut $house_keeping_genes $genome
    """

}