process snp_footprint_clustering{
  debug true
  //errorStrategy 'ignore'
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  //maxRetries 3
  label 'process_medium'

  tag "All Samples"   

  //Docker Image
  container 'prc992/snp_dendrogram:v1.1'
  //publishDir "$path_sample_snp_footprint", mode : 'copy'

  input:
  path(vcfGzFiles)
  path (chRSNPFootprint)
  
  //exec:
  //path_sample_snp_footprint = path_analysis + "/snp_fingerprint" 

  //output:
  //path ('*.pdf')

  //Rscript $chRSNPFootprint

  script:
  """
  Rscript $chRSNPFootprint
  """

}
