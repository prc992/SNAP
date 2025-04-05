process moveSoftFiles {
    label 'low_cpu_low_mem'
    tag "All Samples"
    
    input:
    val(_)
    
    script:
    """
    set -euo pipefail
    mkdir -p ${workflow.projectDir}/${params.outputFolder}/software_versions || true
    mkdir -p ${workflow.projectDir}/${params.outputFolder}/stats_files || true

    echo "Moving mqc_versions.yml files to ${workflow.projectDir}/${params.outputFolder}/software_versions"
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*mqc_versions.yml' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/software_versions/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*AfterFilter*' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*lc_extrap*' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*flagstat' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*idxstats' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*stats' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*metrics.txt' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder}/reports/multiqc/  -type f ! -name "*.html" -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*gz_trimming_report.txt' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name '*-dummy.txt' -exec rm {} \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name 'mark_for_deletion_*' -exec rm {} \\; || true
    find ${workflow.projectDir}/${params.outputFolder} -type f -name 'igv_housekeeping_genes_mqc.html' -exec mv {} ${workflow.projectDir}/${params.outputFolder}/stats_files/ \\; || true

    echo "Finalizing moveSoftFiles step" >&1
    echo "Stub output for testing" > .command.out
    """
}