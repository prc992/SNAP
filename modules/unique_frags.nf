process unique_frags {
  label 'process_low'

  container = params.containers.ubuntu

  tag "Sample - $sampleId" 
  publishDir "$path_sample_peaks", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path (sampleBed),val(_)

  output:
  path ('*.csv')

  exec:
  path_sample_peaks = path_analysis + "/peaks/" + sampleId
  strCSV = sampleId + '_unique_frags.csv'
  
  script:
  """
  echo $sampleId && wc -l $sampleBed | cut -f1 -d' '  >> $strCSV
  """

}
