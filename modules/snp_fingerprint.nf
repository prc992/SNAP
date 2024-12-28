process snp_fingerprint {
  label 'process_medium'
  
  //Docker Image
  container = 'quay.io/mcrotti1/bcftools'

  tag "Sample - $sampleId" 
  publishDir "$path_sample_snp_fingerprint", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sampleBam),path (sampleBai),val(_)
  each path (snps_ref)
  each path (file_fa)


  exec:
  path_sample_snp_fingerprint = path_analysis + "/snp_fingerprint/" + sampleId
  strVCFgz = sampleId + '.vcf.gz'

  output:
  path("*.vcf.gz")

  script:
  """
  bcftools mpileup --threads $task.cpus -Ou -R $snps_ref -f $file_fa $sampleBam | bcftools call --threads $task.cpus -c | bgzip --threads $task.cpus > $strVCFgz
  """

}
