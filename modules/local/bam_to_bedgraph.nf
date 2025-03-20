process bam_to_bedgraph {
  label 'med_cpu_med_mem'
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId),val(control),path(sampleBam),path (indexBam),path ("bam_to_bedgraph_mqc_versions.yml")

  exec:
  strbedgraph = sampleId + '.bedgraph'

  output:
  tuple val(sampleId),val(control),path('*.bedgraph'),path ("bam_to_bedgraph_mqc_versions.yml")

  script:
  """
  bedtools genomecov -ibam $sampleBam -bg > $strbedgraph

  cat <<-END_VERSIONS > bam_to_bedgraph_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}