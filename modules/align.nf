process align {
  label 'high_cpu_high_mem'

  // Docker Image
  container = 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  tag "Sample - $sampleId"
  publishDir "$path_sample_align", mode: 'copy'

  input:
  tuple path(file1), path(file2)
  tuple val(sampleId), val(path), path(_), path(_)
  each path (genomeFile)
  each path (genomeIndexFiles)

  output:
  path("*.bam")

  exec:
  String strBam = sampleId + '.bam'
  path_sample_align = path + "/align/" + sampleId

  script:
  """
  # Find the reference genome file from the input
  file_fa=\$(find -L . -type f -name "*.fa")
  if [ -z "\$file_fa" ]; then
    echo "Error: No .fa file found in the provided path"
    exit 1
  fi
  echo "Reference genome found: \$file_fa"

  # Find the reference index file from the input
  file_pac=\$(find -L . -type f -name "*.fa.pac")
  if [ -z "\$file_fa" ]; then
    echo "Error: Index not found. Please create index for reference genome."
    exit 1
  fi
  echo "Index file found: \$file_pac"

  # Perform alignment
  echo "Running alignment with bwa mem..."
  bwa mem $genomeFile $file1 $file2 -t $task.cpus | \
  samtools view --threads $task.cpus -Sb -u > $strBam
  """
}