process fastqc {
  label 'med_cpu_med_mem'

  container = params.containers.fastqc

  tag "Sample - $sampleId"  
  publishDir "${workflow.projectDir}/${params.outputFolder}/fastqc/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId), val(enrichment_mark),path(read1), path(read2)

  output:
  tuple val(sampleId), path ('*_fastqc.html'),path ('*_fastqc.zip'),path ("fastqc_mqc_versions.yml")
  
  script:
  """
  if [ -z "$read2" ]; then
      # Single-end
      fastqc $read1 --threads $task.cpus
  else
      # Paired-end
      fastqc --threads $task.cpus $read1 $read2
  fi

  cat <<-END_VERSIONS > fastqc_mqc_versions.yml
  "${task.process}":
     fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
  END_VERSIONS
  """


}
