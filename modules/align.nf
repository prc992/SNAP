process align {
  debug true
  label 'high_cpu_high_mem'

  // Docker Image
  container = 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  tag "Sample - $sampleId"
  //publishDir "$path_sample_align", mode: 'copy'

  input:
  tuple val(sampleId), path(trimmedFile) from chTrimFiles

  /*input:
  path(trimmed_files)
  each path (genomeFile)
  each path (genomeIndexFiles)*/

  //exec:
  //String strBam = sampleId + '.bam'
  //path_sample_align = path + "/align/" + sampleId

  script:
  """
  # Debugging: Print input paths
  echo "Trimmed files: $trimmed_files"
  echo "sampleId : $sampleId"
  
  # Print number of files in trimmed_files
  num_files=\$(ls -1 ${trimmed_files} | wc -l)
  echo "Number of files in trimmed_files: \$num_files"
  """
}