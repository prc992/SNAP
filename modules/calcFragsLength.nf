process calcFragsLength {
  label 'med_cpu_med_mem'
  container = params.containers.deeptools

  tag "Sample - $sampleId"  
  publishDir "$path_sample_frags", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sortedBam),path (sampleBamIndex),val(_)

  output:
  tuple path("*fragment_sizes.txt"), path ("bamPEFragmentSize_mqc_versions.yml")

  exec:
  path_sample_frags = path_analysis + "/frag/" + sampleId

  script:
  """
  bamPEFragmentSize -b $sortedBam --outRawFragmentLengths ${sampleId}.fragment_sizes.txt

  cat <<-END_VERSIONS > bamPEFragmentSize_mqc_versions.yml
    "${task.process}":
    bamPEFragmentSize: \$(bamPEFragmentSize --version | sed -e "s/bamPEFragmentSize //g")
  END_VERSIONS
  """
}