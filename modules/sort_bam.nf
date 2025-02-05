process sort_bam {
  label 'low_cpu_high_mem'

  container = params.containers.samtools

  tag "Sample - $sampleId" 
  publishDir "${workflow.projectDir}/${params.outputFolder}/align/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),path(sampleBam),val(_)
  
  output:
  tuple val(sampleId),path('*.bam'),path ("samtools_sort_mqc_versions.yml")

  exec:
  String strBam = sampleId + '.sorted.bam'

  script:
  """
  samtools sort -@ $task.cpus $sampleBam -o $strBam

  cat <<-END_VERSIONS > samtools_sort_mqc_versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
