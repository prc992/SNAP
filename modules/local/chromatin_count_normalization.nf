process chromatin_count_normalization {
  label 'med_cpu_med_mem'
  container params.containers.chromatin_count_normalization

  tag "Sample - $sampleId"
  publishDir "${workflow.projectDir}/${params.outputFolder}/chromatin_count_normalization/${sampleId}", mode: 'copy'

  input:
  tuple val(sampleId), val(enrichment_mark), val(read_method), val(_), val(_), path(narrowPeakFile), val(_), val(_)
  tuple val(_), val(_), val(_), val(_), path(bedFile), val(_)
  path (referenceSitesFile)

  output:
  path "output", type: 'dir' 

  script:
  // If referenceSitesFile exist, use as an argument
  def ref_arg = referenceSitesFile ? "--reference-sites ${referenceSitesFile}" : ""

  """
  Rscript /workspace/chromatin_count_norm_v2.R \
    --sample-name ${sampleId} \
    --fragment-file ${bedFile} \
    --target-sites ${narrowPeakFile} \
    ${ref_arg} \
    --verbose
  """
}