process bedGraphToBigWig {
  label 'process_medium'

  //Docker Image
  container = "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path (treat_pileup_bdg),path (control_lambda_bdg),path (_)
  each path (RefGenSizes)

  exec:
  path_sample_peaks = path_analysis + "/peaks/" + sampleId
  bdgFile1 = treat_pileup_bdg
  bdgFile2 = control_lambda_bdg
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
