process trim {
  label 'med_cpu_high_mem'

  container = params.containers.trim_galore

  tag "Sample - $sampleId"
  //publishDir "${workflow.projectDir}/${params.outputFolder}/trim/${sampleId}", mode: 'copy'

  input:
  tuple val(sampleId), val(enrichment_mark), val(control), val(read_method), path(reads)

  output:
  tuple val(sampleId), val(enrichment_mark), val(control), val(read_method), path('*.fq.gz'), path("*report.txt"), path("trim_mqc_versions.yml")

  script:
  def trimCommand = ""

  if (read_method == "PE") {
    trimCommand = "trim_galore --paired ${reads[0]} ${reads[1]} --gzip --cores $task.cpus $params.trimming_params"
  } else {
    trimCommand = "trim_galore ${reads[0]} --gzip --cores $task.cpus $params.trimming_params"
  }

  """
  echo "Running Trim Galore for sample $sampleId in $read_method mode"
  echo "Command executed: $trimCommand"
  $trimCommand

  cat <<-END_VERSIONS > trim_mqc_versions.yml
  "${task.process}":
      trimgalore: \$(trim_galore --version 2>&1 | sed 's/^.*version //; s/Last.*\$//')
      cutadapt: \$(cutadapt --version)
  END_VERSIONS
  """
}