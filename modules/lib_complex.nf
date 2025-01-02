process lib_complex {
  label 'med_cpu_high_mem'
  container = params.containers.picard

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sortedBam),val(_)

  output:
  tuple val(sampleId),path("*.metrics.txt"),path ("picard_EstimateLibraryComplexity_mqc_versions.yml")

  exec:
  String strLib = sampleId + '.LibComplexity.metrics.txt'
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  picard EstimateLibraryComplexity I=$sortedBam  O=$strLib

  cat <<-END_VERSIONS > picard_EstimateLibraryComplexity_mqc_versions.yml
    "${task.process}":
        picard: \$(echo \$(picard EstimateLibraryComplexity --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
  END_VERSIONS
  """
}
