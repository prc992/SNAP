process call_peaks{
  label 'low_cpu_low_mem'

  container = params.containers.macs2

  tag "Sample - $sampleId"   
  
  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId),path(sampleBam),val(_)

  output:
  tuple val(sampleId),path ('*treat_pileup.bdg'),path ('*control_lambda.bdg'),path ('*narrowPeak'),path("*.xls"),path("macs2_mqc_versions.yml")
  
  script:
  """
  macs2 \\
  callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel \\
  -t $sampleBam \\
  -n $sampleId --bdg

  cat <<-END_VERSIONS > macs2_mqc_versions.yml
    "${task.process}":
        macs2: \$(macs2 --version | sed -e "s/macs2 //g")
  END_VERSIONS
  """
}
