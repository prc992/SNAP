/*process multiqc {
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
}*/

process multiqc_fragments_processing {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
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

process multiqc_bam_processing {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
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

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
    input:
    val(chInitReport)
    path (chFilesReport)
    path (chMultiQCConfig)

    output:
    path ("*.html")

    script:
    """
    multiqc . 
    """
}

process multiqc_initialization {
    label 'process_medium'
    tag "All Samples" 
    container = params.containers.multiqc 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
    input:
    path (chFilesReport)

    output:
    path ("*.*")

    script:
    """
    multiqc . 
    """
}