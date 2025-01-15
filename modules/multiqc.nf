process multiqc {
    label 'med_cpu_high_mem'
    container = params.containers.multiqc 
    publishDir "$path_sample_multiqc", mode : 'copy'
    tag "All Samples" 
    
    input:
    tuple val(_),val(_),val (_),val (_),path (bedGraphToBigWig_mqc_versions)
    path (chIGVReport)
    val(_)
    path (chFragmentsSizeFiles)
    //tuple path(chFragmentsSize), val (_)
    path (chFootPrintPDF)
    path (chFragAndPeaks)
    path (chEnrichmentFiles)
    path (configFile)
    //path (chMultiQCHousekeepingReport)
    tuple val(sampleId), val(enrichment_mark),path(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    path ("multiqc_report.html")
    path ("multiqc_data/*.*")
    path ("multiqc_plots/*.*")

    script:
    """
    multiqc . 
    """
}