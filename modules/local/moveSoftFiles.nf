process moveSoftFiles {
    label 'low_cpu_low_mem'
    tag "All Samples"
    
    input:
    val(_)
    
    script:
    """
    mkdir -p ${workflow.projectDir}/${params.outputFolder}/software_versions
    mkdir -p ${workflow.projectDir}/${params.outputFolder}/stats_files

    echo "Moving mqc_versions.yml files to ${workflow.projectDir}/${params.outputFolder}/software_versions"
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*mqc_versions.yml' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/software_versions/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*AfterFilter*' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*lc_extrap*' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*flagstat' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*idxstats' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*stats' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*metrics.txt' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder}/reports/multiqc/  -type f ! -name "*.html" -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*gz_trimming_report.txt' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*-dummy.txt' -exec rm {} \\;
    find ${workflow.projectDir}/${params.outputFolder} -type f -name 'igv_housekeeping_genes_mqc.html' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\;
    """
}