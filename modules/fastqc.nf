process fastqc {
  label 'med_cpu_med_mem'
    
  //Docker Image
  container = 'quay.io/biocontainers/fastqc:0.11.9--0'

  tag "Sample - $sampleId"  
  publishDir "$path_sample_fastqc", mode : 'copy'
  
  input:
  tuple val(sampleId), val(enrichment_mark),val(path),path(read1), path(read2)

  exec:
  path_sample_fastqc = path + "/fastqc/" + sampleId

  output:
  path ('*_fastqc.html')
  path ('*_fastqc.zip')
  tuple val(sampleId),path ("fastqc_mqc_versions.yml")
  
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
