process bam_to_bed {
  label 'process_medium'

  //Docker Image
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),val(control),path(sampleBam),val(_),val(_)

  exec:
  String strBed = sampleId + '.bed'

  output:
  tuple val(sampleId),val(control),path ('*.bed'),path ("bam_to_bed_mqc_versions.yml")

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
