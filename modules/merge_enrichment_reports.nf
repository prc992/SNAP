process merge_enrichment_reports {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "Sample - $sampleId"

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/enrichment/", mode: 'copy'

    input:
    path (chEnrichmentFilesReport)
    each path (chMultiQCEnrichmentHeader)
    each path (chMergeReportEnrichment)
    tuple val(sampleId), val(enrichment_mark), val(read1), val(read2)

    output:
    path ("*.*")

    script:
    """
    python $chMergeReportEnrichment
    
    # If no file *_mqc.csv was created create an empty one
    if ! ls *_mqc.csv 1> /dev/null 2>&1; then
        touch ${sampleId}.csv
    fi
    """
}