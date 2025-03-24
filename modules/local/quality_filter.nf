process quality_filter {
    label 'low_cpu_low_mem'
    container = params.containers.samtools
    
    tag "Sample - $sampleId" 
    publishDir "${workflow.projectDir}/${params.outputFolder}/align/${sampleId}", mode : 'copy'

    input:
    tuple val(sampleId),val(control),path(sampleBam),val(_)
    tuple val(sampleId), val(_),val(_),val(reads)

    exec:
    String strBam = sampleId + '.filtered.unique.sorted.bam'
  
    output:
    tuple val(sampleId),val(control),path('*.bam'),path ("samtools_QualityFilter_mqc_versions.yml")

    script:
    def singlePairedFilter = reads.size() > 1 ? "-f $params.filter_samtools.inclusion_flag" : ""

    """
    samtools view -bh $singlePairedFilter -F $params.filter_samtools.exclusion_flag -q $params.filter_samtools.min_qc --threads $task.cpus $sampleBam > $strBam

    cat <<-END_VERSIONS > samtools_QualityFilter_mqc_versions.yml
    "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}