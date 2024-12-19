process fetch_chrom_sizes{
  label 'process_low'
  
  //Docker Image
  container = 'quay.io/biocontainers/ucsc-fetchchromsizes:377--ha8a8165_3'

  tag "Fetch Chrom Sizes - $genome"   
  //publishDir "$refDir/genome", mode : 'copy'

  input:
  val genome
  path refDir

  output:
  path ('*.sizes')
  
  script:
  def chromSizeFile = "${genome}.chrom.sizes"
  """
  if [ ! -f ${refDir}/${chromSizeFile} ]; then
    fetchChromSizes $genome > ${chromSizeFile}
  else
    echo "File ${refDir}/${chromSizeFile} already exists. Skipping download."
  fi
  
  ln -s ${refDir}/${chromSizeFile} ${chromSizeFile}
  
  """
}
