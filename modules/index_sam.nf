process index_sam {
  label 'low_cpu_low_mem'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "Sample - $sampleId"   
  publishDir "$path_sample_align", mode : 'copy'
  
  input:
  tuple path (sampleBam), path(_)
  tuple val (sampleId), val(path), path(_), path(_)

  exec:
  path_sample_align = path + "/align/" + sampleId

  output:
  path ('*.bai')

  script:
  """
  samtools index -@ $task.cpus $sampleBam
  """
}
