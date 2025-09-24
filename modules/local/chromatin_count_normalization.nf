process chromatin_count_normalization_batch {
  label 'med_cpu_med_mem'
  container params.containers.chromatin_count_normalization

  tag "Sample - $sampleId"
  publishDir "${workflow.projectDir}/${params.outputFolder}/chromatin_count_normalization/batch", mode: 'copy'

  input:
  path (bedFiles)
  each path (referenceSitesFile)
  each path (targetSitesFile)

  output:
  path "output", type: 'dir' 

  script:
  // If referenceSitesFile exist, use as an argument
  def ref_arg = params.chromatin_count_reference ? "--reference-sites ${referenceSitesFile}" : ""

  """
  echo "Total de BEDs: ${bedFiles.size()}"
  for f in ${bedFiles}; do
    echo \$f
  done

  echo "Gerando arquivo sample_name com header..."
  > sample_name
  echo "sample_name" >> sample_name        # cabeçalho
  for f in ${bedFiles}; do
    bn=\$(basename "\$f")                  # ex: DV_53_M2413... .bed
    name="\${bn%.bed}"                     # remove extensão .bed
    printf "%s\\n" "\$name" >> sample_name
  done

  echo "Arquivo sample_name:"
  cat sample_name

  Rscript /workspace/chromatin_count_norm_v2.R \
  --samplesheet sample_name \
  --target-sites ${targetSitesFile} \
  --frags-dir . \
  ${ref_arg} \
  --verbose
  """
}

process chromatin_count_normalization_single {
  label 'med_cpu_med_mem'
  container params.containers.chromatin_count_normalization

  tag "Sample - $sampleId"
  publishDir "${workflow.projectDir}/${params.outputFolder}/chromatin_count_normalization/${sampleId}", mode: 'copy'

  input:
  tuple val(sampleId), val(_), val(_), val(_), path(bedFile), val(_)
  path (referenceSitesFile)
  path (targetSitesFile)

  output:
  path "output", type: 'dir' 

  script:
  // If referenceSitesFile exist, use as an argument
  def ref_arg = params.chromatin_count_reference ? "--reference-sites ${referenceSitesFile}" : ""

  """
  Rscript /workspace/chromatin_count_norm_v2.R \
    --sample-name ${sampleId} \
    --fragment-file ${bedFile} \
    --target-sites ${targetSitesFile} \
    ${ref_arg} \
    --verbose
  """
}