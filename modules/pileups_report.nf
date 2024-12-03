process pileups_report{
  errorStrategy 'ignore'
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  //maxRetries 3
  label 'process_medium'

  tag "Sample - $sampleId"   

  //Docker Image
  container ='prc992/pileups-report:v1.1'
  publishDir "$path_sample_pile_ups", mode : 'copy'

  input:
  tuple val(sampleId), val(path),path(read1), path(read2)
  path (chChromSizes)
  tuple path (control_pileup_bw),path (treat_pileup_bw)
  each path (chBED)
  each path (chRPileups)

  exec:
  //Here we extract the reference genome from the param align_ref
  //using string manipulation
  path_sample_pile_ups = path + "/pile_ups/" + sampleId

  //strAlign = '"$params.align_ref"'
  //fim = strAlign.lastIndexOf('/')
  //refGenome = strAlign.substring(fim-4,fim)

  output:
  path ('*.pdf')

  script:
  """
  Rscript $chRPileups $treat_pileup_bw $chBED $chChromSizes $params.genome_ref
  """
}
