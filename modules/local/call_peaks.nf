process call_peaks {
  label 'low_cpu_low_mem'
  container = params.containers.macs2

  tag "Sample - $sampleId"   

  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId), val(enrichment_mark), val(read_method), path(sampleBam), path(sampleControl)

  output:
  tuple val(sampleId), val(enrichment_mark), val(read_method), path('*treat_pileup.bdg'), path('*control_lambda.bdg'), path('*narrowPeak'), path("*.xls"), path("macs2_mqc_versions.yml")

  script:
  def bamFormat = ""
  if (read_method == "PE") {
    bamFormat = "BAMPE"
  } else {
    bamFormat = "BAM"
  }

  """
  echo "Running MACS2 callpeak for sample $sampleId in $read_method mode (format: $bamFormat)"

  if [ ! -s ${sampleControl} ]; then
    macs2 \\
    callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -f $bamFormat --extsize 146 --nomodel \\
    -t $sampleBam \\
    -n $sampleId --bdg
  else
    macs2 \\
    callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -f $bamFormat --extsize 146 --nomodel \\
    -t $sampleBam \\
    -c $sampleControl \\
    -n $sampleId --bdg
  fi

  cat <<-END_VERSIONS > macs2_mqc_versions.yml
  "${task.process}":
      macs2: \$(macs2 --version | sed -e "s/macs2 //g")
  END_VERSIONS
  """
}