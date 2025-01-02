process bam_to_bed {
  label 'process_medium'

  //Docker Image
  container = params.containers.bedtools

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'
  
  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam),val(_)

  exec:
  String strBed = sampleId + '.bed'
  path_sample_peaks = path_analysis + "/peaks/" + sampleId

  output:
  tuple val(sampleId),val(path_analysis),path ('*.bed'),path ("bam_to_bed_mqc_versions.yml")

  script:
  """
  bedtools bamtobed -i \\
  $sampleBam -bedpe 2> /dev/null | \\
  awk 'BEGIN{{OFS="\t";FS="\t"}}(\$1==\$4){{print \$1, \$2, \$6}}' > $strBed

  cat <<-END_VERSIONS > bam_to_bed_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}
