process trim {
  label 'med_cpu_high_mem'

  //Docker Image
  container = "quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"

  tag "Sample *.fq - $sampleId"
  publishDir "$path_sample_trim", mode: 'copy'

  input:
  tuple val(sampleId), val(path),path(read1), val(read2)

  exec:
  path_sample_trim = path + "/trim/" + sampleId

  output:
  path('*.fq.gz')

  script:
  """
  if [ -z "$read2" ]; then
    # Single-end
    trim_galore $read1 --gzip --cores $task.cpus
  else
    # Paired-end
    trim_galore --paired $read1 $read2 --gzip --cores $task.cpus
  fi
  """
}
