process dac_exclusion {
  label 'low_cpu_low_mem'

  container = params.containers.bedtools
  
  tag "Sample - $sampleId"  
  publishDir "${workflow.projectDir}/${params.outputFolder}/align/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId),val(control),path(dedupBam),val(_),val(_)
  each path (sampleDAC)

  exec:
  strBam = sampleId + '.dac_filtered.dedup.unique.sorted.bam'
  strTxt = sampleId + '-dummy.txt'

  output:
  tuple val(sampleId),val(control),path(strBam),path(strTxt),path ("dac_exclusion_mqc_versions.yml")

  script:
  """
  touch $strTxt
  bedtools intersect -v -abam $dedupBam -b $sampleDAC > ${sampleId}.tmp
  mv ${sampleId}.tmp $strBam
  
  cat <<-END_VERSIONS > dac_exclusion_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
	"""
}
