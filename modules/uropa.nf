process json_uropa{
  label 'low_cpu_low_mem'
  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'

  container = "ubuntu:noble-20231221"

  input:
  tuple val(sampleId), val(path),path(_), path(_)

  output:
  path ('cfchip.json')

  exec:
  path_sample_peaks = path + "/peaks/" + sampleId

  script:
  """
  BED_FILE=`find -L ./ -name "*.narrowPeak"`
  echo '{"queries": [' >> cfchip.json
  echo '{"feature":"gene","distance":10000,"filter.attribute" : "gene_type","attribute.value" : "protein_coding","feature.anchor":"start"}],' >> cfchip.json
  echo '"show_attributes":["gene_id", "gene_name","gene_type"],   ' >> cfchip.json
  echo '"priority" : "True",' >> cfchip.json
  echo '"gtf": "gencode.v19.annotation.gtf",' >> cfchip.json
  echo '"bed": "\$BED_FILE"}' >> cfchip.json
  """
}


process uropa {
  debug true
  label 'process_medium'
  //Docker Image
  container = "quay.io/biocontainers/uropa:4.0.3--pyhdfd78af_0"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_peaks", mode : 'copy'
    
  input:
  tuple val(sampleId),val(path_analysis),path (bdgFile), path (narrowPeakFile)
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
  """
}
