process moveSoftFiles {
    label 'low_cpu_low_mem'
    tag "All Samples"
    
    input:
    tuple val(fileHtml), val(files_data), val(files_plots)
    tuple val(_), val(_), val(path_analysis),val(_), val(_)
    
    script:
    """
    mkdir -p ${path_analysis}/software_versions
    echo "Moving mqc_versions.yml files to ${path_analysis}/software_versions"
    find ${path_analysis} -type f -name '*mqc_versions.yml' -exec mv {} ${path_analysis}/software_versions/ \\;
    """
}