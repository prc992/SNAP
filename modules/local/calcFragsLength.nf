process calcFragsLengthDistribuition {
  label 'med_cpu_med_mem'
  container = params.containers.deeptools

  tag "Sample - $sampleId"  
  publishDir "${workflow.projectDir}/${params.outputFolder}/frag/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId),path(sortedBam),path (sampleBamIndex),val(_)

  output:
  tuple path("*fragment_sizes.txt"), path ("bamPEFragmentSize_mqc_versions.yml")

  //bamPEFragmentSize -b $sortedBam --outRawFragmentLengths ${sampleId}.fragment_sizes.txt
  //touch ${sampleId}.fragment_sizes.txt
  script:
  """
  bamPEFragmentSize -b $sortedBam --outRawFragmentLengths ${sampleId}.fragment_sizes.txt

  cat <<-END_VERSIONS > bamPEFragmentSize_mqc_versions.yml
    "${task.process}":
    bamPEFragmentSize: \$(bamPEFragmentSize --version | sed -e "s/bamPEFragmentSize //g")
  END_VERSIONS
  """
}