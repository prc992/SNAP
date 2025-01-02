process align {
  label 'high_cpu_high_mem'

  container = params.containers.bwa_and_samtools

  tag "Sample - $sampleId"
  publishDir "$path_sample_align", mode: 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(trimmedFiles),val(_),val(_)
  each path (genomeFile)
  tuple path (genomeIndexFiles),val(_)
  
  output:
  tuple val(sampleId),val(path_analysis),path('*.bam'),path ("align_mqc_versions.yml")

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
     samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}