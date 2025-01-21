process bam_to_tdf {
  label 'med_cpu_med_mem'
  container = params.containers.igv

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam),path (indexBam),val (_)
  each path(genome)
  each path(genomeIndex)

  exec:
  path_sample_peaks = path_analysis + "/peaks/" + sampleId
  str_tdf = sampleId + '.tdf'

  output:
  tuple val(sampleId),val(path_analysis),path('*.tdf'),path ("bam_to_tdf_mqc_versions.yml")

  script:
  """
  igvtools count $sampleBam $str_tdf $genome

  cat <<-END_VERSIONS > bam_to_tdf_mqc_versions.yml
    "${task.process}":
        igvtools: "\$(igvtools help 2>/dev/null | grep 'IGV Version' | awk '{print \$5}')"
  END_VERSIONS
  """

  
}