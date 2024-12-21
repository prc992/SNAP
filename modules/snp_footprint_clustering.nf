process snp_footprint_clustering{
  label 'process_medium'

  tag "All Samples"   

  //Docker Image
  container 'prc992/snp_dendrogram:v1.3'
  publishDir "$path_sample_snp_footprint", mode : 'copy'

  input:
  path(vcfGzFiles)
  path (chRSNPFootprint)
  val (chOutputDir)

  exec:
  path_sample_snp_footprint =  chOutputDir + "/reports/multiqc/" 

  output:
  path ('*.pdf')

  script:
  """
  Rscript $chRSNPFootprint
  """

}
