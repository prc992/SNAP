process lenght_fragment_dist_step1{
  label 'process_medium'

  //Docker Image
  container ='quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  tag "Sample - $sampleId"  
  publishDir "$path_sample_frag", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam)
  val(_)

  exec:
  path_sample_frag = path_analysis + "/frag/" + sampleId

  output:
  tuple val(sampleId),val(path_analysis),path('*.txt')

  exec:
  strtxt = sampleId + '_fragment_lengths.txt'

  script:
  """
  samtools view --threads $task.cpus $sampleBam | cut -f 9 | awk ' \$1 <= 1000 && \$1 > 0 ' > $strtxt
  """

}

process lenght_fragment_dist_step2{
  label 'process_medium'

  //Docker Image
  container ='pegi3s/r_data-analysis'

  tag "$sampleId" 
  publishDir "$path_sample_frag", mode : 'copy'

  output:
  path ('*.png')

  input:
  tuple val(sampleId),val(path_analysis),path(fragLeng)
  each path (chRfrag_plotFragDist)

  exec:
  strPNG = sampleId + '_fragDist.png' 
  path_sample_frag = path_analysis + "/frag/" + sampleId

  script:
  """
  Rscript $chRfrag_plotFragDist $fragLeng $strPNG $sampleId
  """

}
