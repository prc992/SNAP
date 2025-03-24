process trim {
  label 'med_cpu_high_mem'

  container = params.containers.trim_galore

  tag "Sample - $sampleId"
  publishDir "${workflow.projectDir}/${params.outputFolder}/trim/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId), val(enrichment_mark),val(control),path(reads)
  
  output:
  tuple val(sampleId),val(control),path('*.fq.gz'),path("*report.txt"),path ("trim_mqc_versions.yml")

  script:
  def read1 = reads[0]
  def read2 = reads.size() > 1 ? reads[1] : null

  """
  if [ -n "$read2" ]; then
      echo "Paired-end mode: $read1 + $read2"
      trim_galore --paired $read1 $read2 --gzip --cores $task.cpus
  else
      echo "Single-end mode: $read1"
      trim_galore $read1 --gzip --cores $task.cpus
  fi

  cat <<-END_VERSIONS > trim_mqc_versions.yml
  "\${task.process}":
      trimgalore: \$(trim_galore --version 2>&1 | sed 's/^.*version //; s/Last.*\$//')
      cutadapt: \$(cutadapt --version)
  END_VERSIONS
  """
}
