process signal_report_lite {
    label 'low_cpu_low_mem'
    container = params.containers.python

    tag "All Samples"
    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/metrics_lite/", mode: 'copy'

    input:
    path chReportQualityLite

    output:
    path "QualityMetrics.csv" 

    script:
    """
    touch test.csv
    """
}