process bedGraphToBigWig {
  label 'process_medium'
  container = params.containers.bedgraphtobigwig


  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path (treat_pileup_bdg),path (control_lambda_bdg),path (narrowPeak),path(xlsFile),val(_)
  each path (RefGenSizes)

  exec:
  path_sample_peaks = path_analysis + "/peaks/" + sampleId
  bdgFile1 = treat_pileup_bdg
  bdgFile2 = control_lambda_bdg
  bdgFile1_out = bdgFile1 + ".bw"
  bdgFile2_out = bdgFile2 + ".bw"
  //fileNameOutput = sampleId + "_treat_pileup.bdg.bw"

  output:
  tuple val(sampleId),val(path_analysis),path ("*control_lambda.bdg.bw"),path ("*treat_pileup.bdg.bw"),path ("bedGraphToBigWig_mqc_versions.yml")

  script:
  """
  bedGraphToBigWig $bdgFile1 $RefGenSizes $bdgFile1_out &&
  bedGraphToBigWig $bdgFile2 $RefGenSizes $bdgFile2_out

  cat <<-END_VERSIONS > bedGraphToBigWig_mqc_versions.yml
    "${task.process}":
        ucsc: $params.containers.bedgraphtobigwig_version
  END_VERSIONS
  """

}
