process multiqc {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
    input:
    //path (chIGVReport)
    //path (chFragmentsSizeFiles)
    //path (chFootPrintPDF)
    //path (chPeaks)
    //path (chFrag)
    //path (chEnrichmentFiles)
    path (chFinalTrigger)
    path (chMultiQCConfig)
    path (chAllPreviousFiles)


    output:
    path ("*.*")

    script:
    """
    multiqc . 
    """
}