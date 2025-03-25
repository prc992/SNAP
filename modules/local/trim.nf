process trim {
  label 'med_cpu_high_mem'

  container = params.containers.trim_galore

  tag "Sample - $sampleId"
  publishDir "${workflow.projectDir}/${params.outputFolder}/trim/${sampleId}", mode: 'copy'

  input:
  tuple val(sampleId), val(enrichment_mark), val(control), path(reads)

  output:
  tuple val(sampleId), val(control), val(read_method), path('*.fq.gz'), path("*report.txt"), path("trim_mqc_versions.yml")

  script:
  def read1 = reads[0]
  def read2 = reads.size() > 1 ? reads[1] : null
  def read_method = reads.size() > 1 ? "PE" : "SE"
  def trimCommand = ""

  if (read_method == "PE") {
    trimCommand = "trim_galore --paired $read1 $read2 --gzip --cores $task.cpus"
  } else {
    trimCommand = "trim_galore $read1 --gzip --cores $task.cpus"
  }

  """
  echo "Running Trim Galore for sample $sampleId in $read_method mode"
  $trimCommand

  cat <<-END_VERSIONS > trim_mqc_versions.yml
  "${task.process}":
      trimgalore: \$(trim_galore --version 2>&1 | sed 's/^.*version //; s/Last.*\$//')
      cutadapt: \$(cutadapt --version)
  END_VERSIONS
  """
}