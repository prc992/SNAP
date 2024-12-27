process dedup {
  label 'med_cpu_med_mem'

  //Docker Image
  container = 'quay.io/biocontainers/picard:2.27.4--hdfd78af_0'

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(uniqueBam)
  val(_)

  exec:
  path_sample_align = path_analysis + "/align/" + sampleId
  strDedupBam = sampleId + '.dedup.unique.sorted.bam'
  strTxt = sampleId + '-MarkDuplicates.metrics.txt'

  output:
  tuple val(sampleId),val(path_analysis),path('*.bam'),path("*.txt")

  script:
  """
  picard MarkDuplicates \\
  I=$uniqueBam \\
  O=$strDedupBam \\
  REMOVE_DUPLICATES=true \\
  ASSUME_SORT_ORDER=coordinate \\
  VALIDATION_STRINGENCY=LENIENT \\
  METRICS_FILE=$strTxt
	"""
}
