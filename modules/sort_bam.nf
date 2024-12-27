process sort_bam {
  label 'low_cpu_high_mem'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "Sample - $sampleId" 
  publishDir "$path_sample_align", mode : 'copy'
  
  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam)
  val(_)
  
  
  output:
  tuple val(sampleId),val(path_analysis),path('*.bam')
  tuple val(sampleId),path ("samtools_sort_mqc_versions.yml")

  exec:
  String strBam = sampleId + '.sorted.bam'
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  samtools sort -@ $task.cpus $sampleBam -o $strBam

  cat <<-END_VERSIONS > samtools_sort_mqc_versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
