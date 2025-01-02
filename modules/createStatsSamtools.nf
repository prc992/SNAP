process createStatsSamtools {
    label 'low_cpu_low_mem'
    container = params.containers.samtools
    publishDir "$path_sample_align", mode : 'copy'
    
    tag "Sample - $sampleId" 

    input:
    tuple val(sampleId),val(path_analysis),path(sampleBam),val(_)

    exec:
    path_sample_align = path_analysis + "/align/" + sampleId

    output:
    tuple val(sampleId),path ('*.stats'),path ('*.idxstats'),path ('*.flagstat'),path ("samtools_stats_mqc_versions.yml")

    script:
    """
    # Generate stats file
    samtools stats $sampleBam > ${sampleId}.notFiltered.stats

    # Generate idxstats file
    samtools idxstats $sampleBam > ${sampleId}.notFiltered.idxstats

    # Generate flagstat file
    samtools flagstat $sampleBam > ${sampleId}.notFiltered.flagstat

    cat <<-END_VERSIONS > samtools_stats_mqc_versions.yml
    "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process createStatsSamtoolsfiltered {
    label 'low_cpu_low_mem'
    container = params.containers.samtools
    publishDir "$path_sample_align", mode : 'copy'
    
    tag "Sample - $sampleId" 

    input:
    tuple val(sampleId),val(path_analysis),path(sampleBam),val(_)

    exec:
    path_sample_align = path_analysis + "/align/" + sampleId

    output:
    tuple val(sampleId),path ('*.stats'),path ('*.idxstats'),path ('*.flagstat'),path ("samtools_stats_filtered_mqc_versions.yml")

    script:
    """
    # Generate stats file
    samtools stats $sampleBam > ${sampleId}.AfterFilter.stats

    # Generate idxstats file
    samtools idxstats $sampleBam > ${sampleId}.AfterFilter.idxstats

    # Generate flagstat file
    samtools flagstat $sampleBam > ${sampleId}.AfterFilter.flagstat

    cat <<-END_VERSIONS > samtools_stats_filtered_mqc_versions.yml
    "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}