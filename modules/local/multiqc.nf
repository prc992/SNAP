process multiqc {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
    input:
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

process multiqc_bam_processing {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/BAM_PROCESSING/", mode : 'copy'
    
    input:
    path (chFilesReport)
    path (chMultiQCConfig)

    output:
    path ("*.*")

    script:
    """
    multiqc . 
    """
}

process multiqc_bam_signal_processing {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/BAM_SIGNAL_PROCESSING/", mode : 'copy'
    
    input:
    path (chFilesReport)
    path (chMultiQCConfig)

    output:
    path ("*.*")

    script:
    """
    multiqc . 
    """
}

process multiqc_initialization {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/INITIALIZATION/", mode : 'copy'
    
    input:
    path (chFilesReport)

    output:
    path ("*.*")

    script:
    """
    multiqc . 
    """
}