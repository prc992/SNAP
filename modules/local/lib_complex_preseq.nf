process lib_complex_preseq {
  label 'med_cpu_high_mem'
  container = params.containers.preseq

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/align/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId),path(sortedBam),val(_)

  output:
  tuple val(sampleId),path("*.lc_extrap.txt"),path ("preseq_mqc_versions.yml")

  script:
  """
  preseq lc_extrap -B $sortedBam > ${sampleId}.lc_extrap.txt

  cat <<-END_VERSIONS > preseq_mqc_versions.yml
  "${task.process}":
    preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
  END_VERSIONS
  """
}