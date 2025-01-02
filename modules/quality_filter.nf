process quality_filter {
    label 'low_cpu_low_mem'
    container = params.containers.samtools
    publishDir "$path_sample_align", mode : 'copy'
    
    tag "Sample - $sampleId" 

    input:
    tuple val(sampleId),val(path_analysis),path(sampleBam),val(_)

    exec:
    String strBam = sampleId + '.filtered.unique.sorted.bam'
    path_sample_align = path_analysis + "/align/" + sampleId

    output:
    tuple val(sampleId),val(path_analysis),path('*.bam'),path ("samtools_QualityFilter_mqc_versions.yml")

    script:
    """
    samtools view -bh -f 3 -F 3844 -q 30 --threads $task.cpus $sampleBam > $strBam

    cat <<-END_VERSIONS > samtools_QualityFilter_mqc_versions.yml
    "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}