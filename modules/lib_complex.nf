process lib_complex {
  label 'med_cpu_high_mem'

  //Docker Image
  container = "quay.io/biocontainers/picard:2.27.4--hdfd78af_0"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  path sampleBam
  tuple val(sampleId), val(path),path(_), path(_)

  output:
  path("*.csv")

  exec:
  String strBam = sampleId + '_lib_complexity.csv'
  path_sample_align = path + "/align/" + sampleId

  script:
  """
  picard EstimateLibraryComplexity I=$sampleBam  O=$strBam
  """
}
