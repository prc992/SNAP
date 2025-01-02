process bam_to_bedgraph {
  label 'med_cpu_med_mem'
  container = params.containers.bedtools

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam),path (indexBam),path ("bam_to_bedgraph_mqc_versions.yml")

  exec:
  path_sample_peaks = path_analysis + "/peaks/" + sampleId
  strbedgraph = sampleId + '.bedgraph'

  output:
  tuple val(sampleId),val(path_analysis),path('*.bedgraph'),path ("bam_to_bedgraph_mqc_versions.yml")

  script:
  """
  bedtools genomecov -ibam $sampleBam -bg > $strbedgraph

  cat <<-END_VERSIONS > bam_to_bedgraph_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}