process bedGraphToBigWig {
  label 'process_medium'

  //Docker Image
  container = "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  path (RefGenSizes)
  tuple val (_), path (bdgFiles)
  tuple val(sampleId), val(path),path(_), path(_)

  exec:
  path_sample_peaks = path + "/peaks/" + sampleId
  bdgFile1 = bdgFiles.first()
  bdgFile2 = bdgFiles.last()
  bdgFile1_out = bdgFile1 + ".bw"
  bdgFile2_out = bdgFile2 + ".bw"
  //fileNameOutput = sampleId + "_treat_pileup.bdg.bw"

  output:
  tuple path ("*control_lambda.bdg.bw"),path ("*treat_pileup.bdg.bw")

  script:
  """
  bedGraphToBigWig $bdgFile1 $RefGenSizes $bdgFile1_out &&
  bedGraphToBigWig $bdgFile2 $RefGenSizes $bdgFile2_out
  """

}
