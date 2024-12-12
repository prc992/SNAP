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
  tuple val(sampleId),val(path_analysis),path (control_lambda_bdg_bw),path (treat_pileup_bdg_bw)
  each path (chChromSizes)
  each path (chBED)
  each path (chRPileups)
  each path (genomeFile)

  exec:
  path_sample_pile_ups = path_analysis + "/pile_ups/" + sampleId

  output:
  path ('*.pdf')

  script:
  """
  Rscript $chRPileups $treat_pileup_bdg_bw $chBED $chChromSizes $genomeFile
  """
}
