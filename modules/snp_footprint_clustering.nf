process snp_footprint_clustering{
  label 'process_medium'

  tag "All Samples"   

  //Docker Image
  container 'prc992/snp_dendrogram:v1.3'
  publishDir "$path_sample_snp_footprint", mode : 'copy'

  input:
  tuple val(sampleId),path(vcfGzFiles),path (yaml)
  path (chRSNPFootprint)
  tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

  exec:
  path_sample_snp_footprint =  path_analysis + "/reports/multiqc/" 

  output:
  path ('*.jpg')

  script:
  """
  Rscript $chRSNPFootprint
  """

}
