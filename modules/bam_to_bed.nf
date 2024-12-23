process bam_to_bed {
  label 'process_medium'

  //Docker Image
  container ='quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'
  
  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam)

  exec:
  String strBed = sampleId + '.bed'
  path_sample_peaks = path_analysis + "/peaks/" + sampleId

  output:
  tuple val(sampleId),val(path_analysis),path ('*.bed')

  script:
  """
  bedtools bamtobed -i \\
  $sampleBam -bedpe 2> /dev/null | \\
  awk 'BEGIN{{OFS="\t";FS="\t"}}(\$1==\$4){{print \$1, \$2, \$6}}' > $strBed
  """
}
