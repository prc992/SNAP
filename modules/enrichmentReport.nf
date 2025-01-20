process enrichmentReport {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "Sample - $sampleId" 

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
    # Print the value of path_sample_multiqc process enrichmentReport 
    echo "The path for MultiQC reports is: ${path_sample_multiqc}"

    python $chReportEnrichment --mark ${enrichment_mark} --samplename ${sampleId}
    """
}