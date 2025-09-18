process chromatin_count_normalization {
  label 'med_cpu_med_mem'
  container params.containers.chromatin_count_normalization
  containerOptions '--platform=linux/arm64'

  tag "Sample - $sampleId"
  publishDir "${workflow.projectDir}/${params.outputFolder}/chromatin_count_normalization/${sampleId}", mode: 'copy'

  input:
  tuple val(sampleId), val(enrichment_mark), val(read_method), val(_), val(_), path(narrowPeakFile), val(_), val(_)
  tuple val(_), val(_), val(_), val(_), path(bedFile), val(_)

  output:
  path "logs",     type: 'dir' 
  path "matrices", type: 'dir' 

  script:
  """
  Rscript /workspace/chromatin_count_norm_v2.R --sample-name ${sampleId} --fragment-file ${bedFile} --target-sites ${narrowPeakFile} --verbose
  """
}