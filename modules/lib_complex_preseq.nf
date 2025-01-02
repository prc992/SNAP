process lib_complex_preseq {
  label 'med_cpu_high_mem'
  container = params.containers.preseq

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sortedBam),val(_)

  output:
  tuple val(sampleId),path("*.lc_extrap.txt"),path ("preseq_mqc_versions.yml")

  exec:
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  preseq lc_extrap -B $sortedBam > ${sampleId}.lc_extrap.txt

  cat <<-END_VERSIONS > preseq_mqc_versions.yml
  "${task.process}":
    preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
  END_VERSIONS
  """
}