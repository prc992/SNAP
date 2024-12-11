process fetch_chrom_sizes{
  label 'process_low'
  
  //Docker Image
  container = 'quay.io/biocontainers/ucsc-fetchchromsizes:377--ha8a8165_3'

  tag "Fetch Chrom Sizes - $genome"   
  //publishDir "$refDir/genome", mode : 'copy'

  input:
  val genome
  path refDir

  exec:
  refGenomeFile = genome + '.chrom.sizes'

  output:
  path ('*.sizes')
  
  script:
  """
  fetchChromSizes $genome > $refGenomeFile
  ln -s ${refDir}/${refGenomeFile} ${refGenomeFile}
  """
}
