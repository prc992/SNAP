process bam_to_bed {
  label 'process_medium'

  //Docker Image
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),val(control),path(sampleBam),val(_),val(_)
  tuple val(sampleId), val(_),val(_),val(reads)

  exec:
  String strBed = sampleId + '.bed'
  def is_paired = reads.size() > 1 ? true : false

  output:
  tuple val(sampleId),val(control),path ('*.bed'),path ("bam_to_bed_mqc_versions.yml")

  script:
  def bedCommand = is_paired ?
    "bedtools bamtobed -i $sampleBam -bedpe | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} (\$1==\$4){print \$1, \$2, \$6}' > $strBed" :
    "bedtools bamtobed -i $sampleBam | awk 'BEGIN{OFS=\"\\t\"} {print \$1, \$2, \$3}' > $strBed"
  """
  $bedCommand

  cat <<-END_VERSIONS > bam_to_bed_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}
