process frags_and_peaks {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "All Samples"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    path (narrowPeakFiles)
    path (chPeakAllFiles)
    each path (chMultiQCFragPeaksHeader)
    each path (chReportFragHist)
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    path ("frags_and_peaks_mqc.csv")

    script:
    """
    python $chReportFragHist
    """
}