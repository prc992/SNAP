process snp_fingerprint {
  label 'process_medium'
  
  //Docker Image
  container = 'quay.io/mcrotti1/bcftools'

  tag "Sample - $sampleId" 
  publishDir "$path_sample_snp_fingerprint", mode : 'copy'

  input:
  tuple path (sampleBam), val (_)
  each path (snps_ref)
  each path (file_fa)
  tuple val(sampleId), val(path),path(_), path(_)
  path (indexFiles)

  exec:
  path_sample_snp_fingerprint = path + "/snp_fingerprint/" + sampleId
  strVCFgz = sampleId + '.vcf.gz'
  
  output:
  tuple val(sampleId), path("*.vcf.gz")

  script:
  """
  FASTA=`find -L ./ -name "*.fa"`
  bcftools mpileup --threads $task.cpus -Ou -R $snps_ref -f \$FASTA $sampleBam | bcftools call --threads $task.cpus -c | bgzip --threads $task.cpus > $strVCFgz
  """

}
