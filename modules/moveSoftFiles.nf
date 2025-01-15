process moveSoftFiles {
    label 'low_cpu_low_mem'
    tag "All Samples"
    
    input:
    val(_)
    val(_)
    val(_)
    tuple val(_), val(_), val(path_analysis),val(_), val(_)
    
    script:
    """
    mkdir -p ${path_analysis}/software_versions
    mkdir -p ${path_analysis}/stats_files

    echo "Moving mqc_versions.yml files to ${path_analysis}/software_versions"
    find ${path_analysis} -type f -name '*mqc_versions.yml' -exec mv {} ${path_analysis}/software_versions/ \\;
    find ${path_analysis} -type f -name '*AfterFilter*' -exec mv {} ${path_analysis}/stats_files/ \\;
    find ${path_analysis} -type f -name '*lc_extrap*' -exec mv {} ${path_analysis}/stats_files/ \\;
    find ${path_analysis} -type f -name '*flagstat' -exec mv {} ${path_analysis}/stats_files/ \\;
    find ${path_analysis} -type f -name '*idxstats' -exec mv {} ${path_analysis}/stats_files/ \\;
    find ${path_analysis} -type f -name '*stats' -exec mv {} ${path_analysis}/stats_files/ \\;
    find ${path_analysis} -type f -name '*gz_trimming_report.txt' -exec mv {} ${path_analysis}/stats_files/ \\;
    find ${path_analysis} -type f -name '*metrics.txt' -exec mv {} ${path_analysis}/stats_files/ \\;
    """
}