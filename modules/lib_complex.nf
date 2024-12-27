process lib_complex {
  label 'med_cpu_high_mem'

  //Docker Image
  container = "quay.io/biocontainers/picard:2.27.4--hdfd78af_0"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sortedBam)
  val(_)

  output:
  path("*.metrics.txt")

  exec:
  String strLib = sampleId + '.LibComplexity.metrics.txt'
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  picard EstimateLibraryComplexity I=$sortedBam  O=$strLib
  """
}
