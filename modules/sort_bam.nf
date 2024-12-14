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

  exec:
  String strBam = sampleId + '.sorted.bam'
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  samtools sort -@ $task.cpus $sampleBam -o $strBam
  """
}
