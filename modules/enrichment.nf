process enrichment {
  label 'process_medium'

  container = params.containers.samtools

  tag "Sample - $sampleId"   
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam),val(_)
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
