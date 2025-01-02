process enrichmentReport {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "All Samples"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(_), val(_)
    path(csvFiles)
    each path (chReportEnrichment)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    path ("*_report.csv")

    script:
    """
    python $chReportEnrichment --mark ${enrichment_mark} --samplename ${sampleId}
    """
}