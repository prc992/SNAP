process merge_enrichment_reports {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "Sample - $sampleId"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    path (chEnrichmentFilesReport)
    each path (chMultiQCEnrichmentHeader)
    each path (chMergeReportEnrichment)
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    path ("*_mqc.csv")

    script:
    """
    python $chMergeReportEnrichment
    """
}