process snp_footprint_clustering{
  //errorStrategy 'ignore'
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  //maxRetries 3
  label 'process_medium'

  tag "All Samples"   

  //Docker Image
  container ='prc992/snp_dendrogram:v1.0'
  publishDir "$path_sample_snp_footprint", mode : 'copy'

  input:
  tuple val(sampleId), val(path_analysis),path(_), path(_)
  path(chRSNPFootprint)
  path (vcfFiles)
  
  exec:
  path_sample_snp_footprint = path_analysis + "/snp_fingerprint" 

  output:
  path ('*.pdf')

  script:
  """
  Rscript $chRSNPFootprint
  """
}
