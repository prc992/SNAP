process uropa {
  label 'process_medium'
  container = params.containers.uropa

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/uropa/${sampleId}", mode : 'copy'
    
  input:
  tuple val(sampleId),path (treat_pileup_bdg),path (control_lambda_bdg),path (narrowPeakFile),path(xlsFile),val(_)
  each path (gtf_file)

  
  output:
  path ('*.*')
  
  script:
  """
  echo '{"queries": [' >> cfchip.json
  echo '{"feature":"gene","distance":10000,"filter.attribute" : "gene_type","attribute.value" : "protein_coding","feature.anchor":"start"}],' >> cfchip.json
  echo '"show_attributes":["gene_id", "gene_name","gene_type"],   ' >> cfchip.json
  echo '"priority" : "True",' >> cfchip.json
  echo '"gtf": "$gtf_file",' >> cfchip.json
  echo '"bed": "$narrowPeakFile"}' >> cfchip.json

  uropa -i cfchip.json -t $task.cpus --summary

  cat <<-END_VERSIONS > uropa_mqc_versions.yml
    "${task.process}":
        uropa: \$(uropa --version)
  END_VERSIONS
  """
}
