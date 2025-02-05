process merge_enrichment_reports {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "Sample - $sampleId"

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'

    input:
    path (chEnrichmentFilesReport)
    each path (chMultiQCEnrichmentHeader)
    each path (chMergeReportEnrichment)
    tuple val(sampleId), val(enrichment_mark),val(read1), val(read2)

    output:
    path ("*_mqc.csv")

    script:
    """
    python $chMergeReportEnrichment
    """
}