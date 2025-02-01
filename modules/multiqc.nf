process multiqc {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 
    publishDir "$path_sample_multiqc", mode : 'copy'
    
    input:
    path (chIGVReport)
    val(_)
    path (chFragmentsSizeFiles)
    path (chFootPrintPDF)
    path (chFragAndPeaks)
    path (chEnrichmentFiles)
    path (configFile)
    tuple val(sampleId), val(enrichment_mark),path(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    path ("*.*")

    script:
    """
    # Print the value of path_sample_multiqc process multiqc
    echo "The path for MultiQC reports is: ${path_sample_multiqc}"

    multiqc . 
    """
}