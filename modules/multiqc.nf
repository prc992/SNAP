process multiqc {
    label 'low_cpu_low_mem'
    container = params.containers.multiqc 
    publishDir "$path_sample_multiqc", mode : 'copy'
    tag "All Samples" 
    
    input:
    tuple val(_),val(_),val (_),val (_),path (bedGraphToBigWig_mqc_versions)
    path (chIGVReport)
    val(_)
    tuple path ("frag_len_hist.txt"),path ("frag_len_mqc.yml")
    path (chFootPrintPDF)
    path (chFragAndPeaks)
    path (chEnrichmentFiles)
    path (configFile)
    path (chMultiQCHousekeepingReport)
    tuple val(sampleId), val(enrichment_mark),path(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    file "multiqc_report.html"
    file "multiqc_data/*"

    script:
    """
    multiqc . 
    find $path_analysis -type f -name '*mqc_versions.yml' -delete
    """
}