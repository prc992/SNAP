process align {
  label 'high_cpu_high_mem'
  
  //Docker Image
  container = 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  tag "Sample - $sampleId" 
  publishDir "$path_sample_align", mode : 'copy'
  
  input:
  tuple path(file1),path(file2)
  tuple val(sampleId), val(path),path(_), path(_)
  each path (file_fa)
  
  output:
  path("*.bam")

  exec:
  String strBam = sampleId + '.bam'
  path_sample_align = path + "/align/" + sampleId

  
  script:
  """
  INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
  bwa mem \$INDEX $file1 $file2 -t $task.cpus | samtools view --threads $task.cpus -Sb -u > $strBam
  """
}
