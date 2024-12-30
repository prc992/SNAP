process uropa {
  label 'process_medium'
  //Docker Image
  container = "quay.io/biocontainers/uropa:4.0.3--pyhdfd78af_0"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'
    
  input:
  tuple val(sampleId),val(path_analysis),path (treat_pileup_bdg),path (control_lambda_bdg),path (narrowPeakFile),path(xlsFile),val(_)
  each path (gtf_file)

  exec:
  path_sample_peaks = path_analysis + "/uropa/" + sampleId
  
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
