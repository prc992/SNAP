process unique_sam {
  label 'process_low'

  container = params.containers.samtools

  tag "Sample - $sampleId" 
  publishDir "$path_sample_align", mode : 'copy'
  
  input:
  tuple val(sampleId),val(path_analysis),path(sortedBam),val(_)

  output:
  tuple val(sampleId),val(path_analysis),path('*.bam'),path ("samtools_unique_mqc_versions.yml")

  exec:
  String strBam = sampleId + '.unique.sorted.bam'
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  samtools view --threads $task.cpus -b -q 1 $sortedBam > $strBam

  cat <<-END_VERSIONS > samtools_unique_mqc_versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
