process trim {
  label 'med_cpu_high_mem'

  //Docker Image
  container = "quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"

  tag "Sample - $sampleId"
  publishDir "$path_sample_trim", mode: 'copy'

  input:
  tuple val(sampleId), val(path),path(read1), path(read2)

  exec:
  path_sample_trim = path + "/trim/" + sampleId

  output:
  tuple path("*1.fq.gz"),path("*2.fq.gz") 

  script:
  """
  trim_galore --paired $read1 $read2 --gzip --cores $task.cpus
  """
}
