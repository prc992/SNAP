process frags_and_peaks {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "Sample - $sampleId" 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    

    input:
    path (narrowPeakFiles)
    path (chPeakAllFiles)
    each path (chMultiQCFragsHeader)
    each path (chMultiQCPeaksHeader)
    each path (chReportFragHist)
    tuple val(sampleId), val(enrichment_mark),val(read1), val(read2)

    output:
    path ("*_mqc.csv")

    script:
    """
    python $chReportFragHist
    """
}