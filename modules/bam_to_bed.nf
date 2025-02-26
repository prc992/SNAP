process bam_to_bed {
  label 'process_medium'

  //Docker Image
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),path(sampleBam),val(_)

  exec:
  String strBed = sampleId + '.bed'
  String strBedPE = sampleId + '.bedpe'


  output:
  tuple val(sampleId),path ('*.bed'),path ('*.bedpe'),path ("bam_to_bed_mqc_versions.yml")

  script:
  """
  bedtools bamtobed -i \\
  $sampleBam -bedpe 2> /dev/null | \\
  awk 'BEGIN{{OFS="\t";FS="\t"}}(\$1==\$4){{print \$1, \$2, \$6}}' > $strBed

  bedtools bamtobed -i \\
  $sampleBam -bedpe 2> /dev/null | \\
  awk 'OFS = "\t" {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $6-$2}' | awk '$11 >=0' > $strBedPE

  cat <<-END_VERSIONS > bam_to_bed_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}
