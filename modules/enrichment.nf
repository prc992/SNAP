process enrichment {
  label 'process_medium'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "Sample - $sampleId"   
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam)
  each path (chEnrichmentScript)


  exec:
  path_sample_peaks = path_analysis + "/peaks/" + sampleId
  strCSV = sampleId + '_enrichment_states.csv'

  output:
  path("*.csv")

  script:
  """
  sh $chEnrichmentScript $sampleBam $params.enrichment_states_ref $sampleId >> $strCSV
  """
}
