process multiqc {
    label 'process_medium'
    tag "All Samples-" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
    input:
    path (chIGVReport)
    path (chFragmentsSizeFiles)
    path (chFootPrintPDF)
    path (chFragAndPeaks)
    path (chEnrichmentFiles)
    path (configFile)
    path (chAllPreviousFiles)


    output:
    path ("*.*")

    script:
    """
    multiqc . 
    """
}