process dac_exclusion {
  label 'low_cpu_low_mem'

  //Docker Image
  container ='quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(dedupBam),val(_)
  val(_)
  each path (sampleDAC)

  exec:
  path_sample_align = path_analysis + "/align/" + sampleId
  strBam = sampleId + '.dac_filtered.dedup.unique.sorted.bam'

  output:
  tuple val(sampleId),val(path_analysis),path('*.bam')

  script:
  """
  bedtools intersect -v -abam $dedupBam -b $sampleDAC > $strBam
	"""
}
