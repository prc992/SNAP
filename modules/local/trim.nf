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
  """
  if [[ \${#reads[@]} -eq 2 ]]; then
      # Paired-end
      trim_galore --paired \${reads[0]} \${reads[1]} --gzip --cores $task.cpus
  else
      # Single-end
      trim_galore \${reads[0]} --gzip --cores $task.cpus
  fi

  cat <<-END_VERSIONS > trim_mqc_versions.yml
  "${task.process}":
      trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
      cutadapt: \$(cutadapt --version)
  END_VERSIONS
  """
}
