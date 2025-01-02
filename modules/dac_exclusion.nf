process dac_exclusion {
  label 'low_cpu_low_mem'

  container = params.containers.bedtools
  
  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(dedupBam),val(_),val(_)
  each path (sampleDAC)

  exec:
  path_sample_align = path_analysis + "/align/" + sampleId
  strBam = sampleId + '.dac_filtered.dedup.unique.sorted.bam'

  output:
  tuple val(sampleId),val(path_analysis),path('*.bam'),path ("dac_exclusion_mqc_versions.yml")

  script:
  """
  bedtools intersect -v -abam $dedupBam -b $sampleDAC > $strBam

  cat <<-END_VERSIONS > dac_exclusion_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
	"""
}
