process bam_to_bed {
  label 'process_medium'

  //Docker Image
  container ='quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'
  
  input:
  tuple path (sampleBam), path(_)
  tuple val(sampleId), val(path),path(_), path(_)

  exec:
  String strBed = sampleId + '.bed'
  path_sample_peaks = path + "/peaks/" + sampleId

  output:
  path ('*.bed')

  script:
  """
  bedtools bamtobed -i \\
  $sampleBam -bedpe 2> /dev/null | \\
  awk 'BEGIN{{OFS="\t";FS="\t"}}(\$1==\$4){{print \$1, \$2, \$6}}' > $strBed
  """
}
