process quality_report_lite {
    label 'low_cpu_low_mem'
    container = params.containers.python

    tag "All Samples"
    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/metrics_lite/", mode: 'copy'

    input:
    path chReportQualityLite
    path enrichment_files
    path peaks_files
    path frags_process_report
    path files_fragle_report

    output:
    path "QualityMetrics.csv4" 

    script:
    """
    python $chReportQualityLite 
    """
}