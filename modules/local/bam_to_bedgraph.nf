process bam_to_bedgraph {
  label 'med_cpu_med_mem'
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'

  input:
  tuple val(sampleId),val(enrichment_mark),val(control),val(read_method),path(sortedBam),path (sampleBamIndex),val (_)

  exec:
  strbedgraph = sampleId + '.bedgraph'

  output:
  tuple val(sampleId),val(enrichment_mark),val(control),val(read_method),path('*.bedgraph'),path ("bam_to_bedgraph_mqc_versions.yml")


  script:
  """
  bedtools genomecov -ibam $sortedBam -bg > $strbedgraph

  cat <<-END_VERSIONS > bam_to_bedgraph_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}