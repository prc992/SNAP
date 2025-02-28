process enrichmentReport {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "Sample - $sampleId" 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'

    input:
    tuple val(sampleId), val(enrichment_mark),val(_), val(_)
    path(csvFiles)
    each path (chReportEnrichment)


    output:
    path ("*_report.csv")

    // touch ${sampleId}_report.csv
    script:
    """
    echo "$enrichment_mark"
    python $chReportEnrichment --mark ${enrichment_mark} --samplename ${sampleId}
    """
}