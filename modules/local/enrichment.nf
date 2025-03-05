process enrichment {
  label 'process_medium'

  container = params.containers.samtools

  tag "Sample - $sampleId"   

  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId),path(sampleBam),val(_)
  each path (chEnrichmentScript)


  exec:
  strCSV = sampleId + '_enrichment_states.csv'

  output:
  path("*.csv")

  script:
  """
  sh $chEnrichmentScript $sampleBam $params.enrichment_states_ref $sampleId >> $strCSV
  """
}
