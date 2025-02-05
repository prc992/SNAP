process snp_footprint_clustering{
  label 'process_medium'
  tag "All Samples"   
  container = params.containers.snp_dendrogram

  publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'

  input:
  path(vcfGzFiles)
  path (chRSNPFootprint)
  tuple val(sampleId), val(enrichment_mark),val(read1), val(read2)

  output:
  path ('*.jpg')

  script:
  """
  Rscript $chRSNPFootprint
  """

}
