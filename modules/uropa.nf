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
  //publishDir "$path_sample_peaks", mode : 'copy'
    
  input:
  tuple val(sampleId),val(path_analysis),path (_), path (bedFiles)
  each path (gtf_file)

  exec:
  path_sample_peaks = path_analysis + "/peaks/" + sampleId
  
  //output:
  //path ('*finalhits.bed')
  
  script:
  """
  BED_FILE=`find -L ./ -name "*.narrowPeak"`
  echo '"bed": "\$BED_FILE"}' 
  """
  //uropa -i $json_file -t $task.cpus --summary
}
