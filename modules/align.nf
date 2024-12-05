process align {
  label 'high_cpu_high_mem'

  // Docker Image
  container = 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:8110a70be2bfe7f75a2ea7f2a89cda4cc7732095-0'

  tag "Sample - $sampleId"
  publishDir "$path_sample_align", mode: 'copy'

  input:
  tuple path(file1), path(file2)
  tuple val(sampleId), val(path), path(_), path(_)
  each path (align_ref)

  output:
  path("*.bam")

  exec:
  String strBam = sampleId + '.bam'

  script:
  """
  # Find the reference genome file from the input
  file_fa=$(find -L . -type f -name "*.fa")
  if [ -z "\$file_fa" ]; then
    echo "Error: No .fa file found in the provided path $align_ref"
    exit 1
  fi
  echo "Reference genome found: \$file_fa"

  # Define the index prefix based on the reference file
  INDEX=\$(basename "\$file_fa" .fa)
  echo "Using index prefix: \$INDEX"

  # Check if the index exists, and create it if not
  if [ ! -f "\${INDEX}.bwt" ]; then
    echo "Index not found. Creating index for reference genome: \$file_fa"
    bwa index "\$file_fa"
    if [ \$? -ne 0 ]; then
      echo "Error: Failed to create index for \$file_fa"
      exit 1
    fi
  else
    echo "Index found. Skipping index creation."
  fi

  # Perform alignment
  echo "Running alignment with bwa mem..."
  bwa mem "\$file_fa" "\$file1" "\$file2" -t $task.cpus | \
  samtools view --threads $task.cpus -Sb -u > $strBam

  if [ \$? -eq 0 ]; then
    echo "Alignment completed successfully. Output: $strBam"
  else
    echo "Error: Alignment failed for sample: $sampleId"
    exit 1
  fi
  """
}