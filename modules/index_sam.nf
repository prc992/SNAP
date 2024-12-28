process index_sam {
  label 'low_cpu_low_mem'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "Sample - $sampleId"   
  publishDir "$path_sample_align", mode : 'copy'
  
  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam)
  val(_)

  exec:
  path_sample_align = path_analysis + "/align/" + sampleId

  output:
  tuple val(sampleId),val(path_analysis),path(sampleBam),path ('*.bai'),path ("index_sam_mqc_versions.yml")

  script:
  """
  samtools index -@ $task.cpus $sampleBam

  cat <<-END_VERSIONS > index_sam_mqc_versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
