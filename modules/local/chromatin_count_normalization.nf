process chromatin_count_normalization {
  label 'med_cpu_med_mem'
  container params.containers.chromatin_count_normalization

  tag "Sample - $sampleId"
  publishDir "${workflow.projectDir}/${params.outputFolder}/chromatin_count_normalization/${sampleId}", mode: 'copy'

  input:
  tuple val(sampleId), val(enrichment_mark), val(read_method), val(_), val(_), path(narrowPeakFile), val(_), val(_)
  tuple val(_), val(_), val(_), val(_), path(bedFile), val(_)

output:
  dir "logs"     , emit: logs
  dir "matrices" , emit: matrices

  script:
  """
  echo "Running chromatin-frags-normalization for sample ${sampleId}"

  chromatin-frags-normalization \
    --sample-name ${sampleId} \
    --fragment-file ${bedFile} \
    --target-sites ${narrowPeakFile} \
    --verbose
  """
}