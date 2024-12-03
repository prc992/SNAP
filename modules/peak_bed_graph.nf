process peak_bed_graph{
  label 'low_cpu_low_mem'

  //Docker Image
  container = 'quay.io/biocontainers/macs2:2.2.7.1--py38h4a8c8d9_3'

  tag "Sample - $sampleId"   
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple path (sampleBam), path(_)
  tuple val(sampleId), val(path),path(_), path(_)

  output:
  tuple path ('*.*'),path ('*.bdg')

  exec:
  path_sample_peaks = path + "/peaks/" + sampleId
  
  script:
  """
  macs2 \\
  callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -f BAMPE --extsize 146 --nomodel \\
  -t $sampleBam \\
  -n $sampleId --bdg
  """
}
