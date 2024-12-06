process dac_exclusion {
  label 'low_cpu_low_mem'

  //Docker Image
  container ='quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  path(sampleBam)
  tuple val(sampleId), val(path),path(_), path(_)
  path (sampleDAC)

  exec:
  path_sample_align = path + "/align/" + sampleId
  strBam = sampleId + '.dac_filtered.dedup.unique.sorted.bam'

  output:
  path("*.bam")

  script:
  """
  bedtools intersect -v -abam $sampleBam -b $sampleDAC > $strBam
	"""
}
