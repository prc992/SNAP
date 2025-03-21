process fastqc {
  label 'med_cpu_med_mem'

  container = params.containers.fastqc

  tag "Sample - $sampleId"  
  publishDir "${workflow.projectDir}/${params.outputFolder}/fastqc/${sampleId}", mode : 'move'
  
  input:
  tuple val(sampleId), val(enrichment_mark),path(read1), path(read2), val(control)

  output:
  path('*fastqc.html'),path('fastqc_mqc_versions.yml'),path('*fastqc.zip')
  
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
