process fastqc {
  label 'med_cpu_med_mem'

  container = params.containers.fastqc

  tag "Sample - $sampleId"  
  publishDir "${workflow.projectDir}/${params.outputFolder}/fastqc/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId), val(enrichment_mark),val(control),path(reads)

  output:
  tuple val(sampleId),path('*fastqc.html'),path('fastqc_mqc_versions.yml'),path('*fastqc.zip')
  
  script:
  // Extract read1 and optional read2 from the reads list
  def read1 = reads[0]
  def read2 = reads.size() > 1 ? reads[1] : null
  
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
