process align {
  label 'high_cpu_high_mem'

  // Docker Image
  container = 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  tag "Sample - $sampleId"
  publishDir "$path_sample_align", mode: 'copy'

  input:
  tuple val(sampleId), val(path_analysis),path(trimmedFiles)
  val(_)
  val(_)
  each path (genomeFile)
  each path (genomeIndexFiles)

  output:
  tuple val(sampleId),val(path_analysis),path('*.bam')
  tuple val(sampleId),path ("align_mqc_versions.yml")

  exec:
  String strBam = sampleId + '.bam'
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  bwa mem $genomeFile $trimmedFiles -t $task.cpus | \
   samtools view --threads $task.cpus -Sb -u > $strBam

  cat <<-END_VERSIONS > align_mqc_versions.yml
  "${task.process}":
     bwa: \$( bwa 2>&1 | grep Version | sed -e "s/Version: //g" )
     samtools: \$( samtools --version | sed -e "s/samtools //g" )
  END_VERSIONS
  """
}