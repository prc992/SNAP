process snp_footprint_clustering{
  label 'process_medium'
  tag "All Samples"   
  container = params.containers.snp_dendrogram
  publishDir "$path_sample_snp_footprint", mode : 'copy'

  input:
  path(vcfGzFiles)
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
