process moveSoftFiles {
    label 'low_cpu_low_mem'
    tag "All Samples"
    
    input:
    val(_)
    
    script:
    """
    set -euo pipefail

    SOFT_DIR="${workflow.projectDir}/${params.outputFolder}/software_versions"
    STATS_DIR="${workflow.projectDir}/${params.outputFolder}/stats_files"
    OUT_DIR="${workflow.projectDir}/${params.outputFolder}"

    mkdir -p "\$SOFT_DIR"
    mkdir -p "\$STATS_DIR"

    echo "Moving mqc_versions.yml files to \$SOFT_DIR"
    find "\$OUT_DIR" -type f -name '*mqc_versions.yml' -exec mv -n {} "\$SOFT_DIR/" \; || true

    echo "Moving stat files to \$STATS_DIR"
    find "\$OUT_DIR" -type f -name '*AfterFilter*' -exec mv -n {} "\$STATS_DIR/" \; || true
    find "\$OUT_DIR" -type f -name '*lc_extrap*' -exec mv -n {} "\$STATS_DIR/" \; || true
    find "\$OUT_DIR" -type f -name '*flagstat' -exec mv -n {} "\$STATS_DIR/" \; || true
    find "\$OUT_DIR" -type f -name '*idxstats' -exec mv -n {} "\$STATS_DIR/" \; || true
    find "\$OUT_DIR" -type f -name '*stats' -exec mv -n {} "\$STATS_DIR/" \; || true
    find "\$OUT_DIR" -type f -name '*metrics.txt' -exec mv -n {} "\$STATS_DIR/" \; || true
    find "\$OUT_DIR" -type f -name '*gz_trimming_report.txt' -exec mv -n {} "\$STATS_DIR/" \; || true
    find "\$OUT_DIR" -type f -name 'igv_housekeeping_genes_mqc.html' -exec mv -n {} "\$STATS_DIR/" \; || true

    if [ -d "\$OUT_DIR/reports/multiqc" ]; then
        find "\$OUT_DIR/reports/multiqc" -type f ! -name "*.html" -exec mv -n {} "\$STATS_DIR/" \; || true
    fi

    # Clean up
    find "\$OUT_DIR" -type f -name '*-dummy.txt' -exec rm -f {} \; || true
    find "\$OUT_DIR" -type f -name 'mark_for_deletion_*' -exec rm -f {} \; || true

    # Optional cleanup
    if [ "${params.cleanup}" == "true" ]; then
        echo "Removing work directory at ${workflow.workDir}"
        rm -rf "${workflow.workDir}" || true
    fi

    echo "Finalizing moveSoftFiles step"
    echo "Stub output for testing" > .command.out
    """
}