process multiqc {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
    input:
    path (chIGVReport)
    val(_)
    path (chFragmentsSizeFiles)
    path (chFootPrintPDF)
    path (chFragAndPeaks)
    path (chEnrichmentFiles)
    path (configFile)
    path (chAllPreviousFiles)
    

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    path ("*.*")

    script:
    """
    multiqc . 
    """
}