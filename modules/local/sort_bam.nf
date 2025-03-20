process sort_bam {
  label 'low_cpu_high_mem'

  container = params.containers.samtools

  tag "Sample - $sampleId" 
  
  //publishDir "${workflow.projectDir}/${params.outputFolder}/align/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),val(control),path(sampleBam),val(_)
  
  output:
  tuple val(sampleId),val(control),path('*.bam'),path ("samtools_sort_mqc_versions.yml")

  exec:
  String strBam = sampleId + '.sorted.bam'

  //samtools sort -@ $task.cpus $sampleBam -o $strBam

  script:
  """
  samtools sort -@ $task.cpus $sampleBam -o $strBam

  cat <<-END_VERSIONS > samtools_sort_mqc_versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}

process sort_readname_bam {
  label 'low_cpu_high_mem'

  container = params.containers.samtools

  tag "Sample - $sampleId" 
  
  input:
  tuple val(sampleId),path(sampleBam),val(_),val(_)
  
  output:
  tuple val(sampleId),path('*.bam'),path ("samtools_sort_mqc_versions.yml")

  exec:
  String strBam = sampleId + '.sorted.bam'

  //samtools sort -@ $task.cpus $sampleBam -o $strBam

  script:
  """
  samtools sort -@ $task.cpus -n $sampleBam -o $strBam

  cat <<-END_VERSIONS > samtools_sort_mqc_versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
