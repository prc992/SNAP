process dedup {
  label 'med_cpu_med_mem'

  container = params.containers.picard

  tag "Sample - $sampleId"  
  
  if (params.intermediate_bams == true) {
    publishDir "${workflow.projectDir}/${params.outputFolder}/align/${sampleId}", mode : 'copy'
  }
  

  input:
  tuple val(sampleId),path(uniqueBam),val(_)

  exec:
  strDedupBam = sampleId + '.dedup.unique.sorted.bam'
  strTxt = sampleId + '-MarkDuplicates.metrics.txt'

  output:
  tuple val(sampleId),path('*.bam'),path("*.txt"),path ("picard_MarkDuplicates_mqc_versions.yml")
  
  script:
  """
  picard MarkDuplicates \\
  I=$uniqueBam \\
  O=$strDedupBam \\
  REMOVE_DUPLICATES=true \\
  ASSUME_SORT_ORDER=coordinate \\
  VALIDATION_STRINGENCY=LENIENT \\
  METRICS_FILE=$strTxt

  cat <<-END_VERSIONS > picard_MarkDuplicates_mqc_versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
  END_VERSIONS
	"""

  
}
