process snp_fingerprint {
  label 'process_medium'
  
  //Docker Image
  container = params.containers.bcftools

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
  tuple val(sampleId),path("*.vcf.gz"),path ("snp_fingerprint_mqc_versions.yml")

  script:
  """
  bcftools mpileup --threads $task.cpus -Ou -R $snps_ref -f $file_fa $sampleBam | bcftools call --threads $task.cpus -c | bgzip --threads $task.cpus > $strVCFgz

  cat <<-END_VERSIONS > snp_fingerprint_mqc_versions.yml
  "${task.process}":
    bcftools: \$(bcftools --version | sed -n '1s/bcftools //p')
  END_VERSIONS
  """

}
