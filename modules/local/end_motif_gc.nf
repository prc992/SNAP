process createMotifGCfile {
  label 'process_medium'

  //Docker Image
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/motifs/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),val(control),path(sampleBam),val(_)
  each path (genomeFile)
  each path (genomeIndexFiles)

  exec:
  String strBed = sampleId + '_frags_gc.bed.bed'
  String strBedPE = sampleId + '.bedpe'
  String strBedFilterPE = sampleId + '_filtered.bedpe'
  String strBPr1 = sampleId + '_' + params.nmer + '_bp_r1.bed'
  String strBPr2 = sampleId + '_' + params.nmer + '_bp_r2.bed'
  String strBPr1FA = sampleId + '_' + params.nmer + 'NMER_bp_r1.fa.bed'
  String strBPr2FA = sampleId + '_' + params.nmer + 'NMER_bp_r2.fa.bed'
  String strBPmotif = sampleId + '_' + params.nmer + 'NMER_bp_motif.bed'
  
  output:
  path ('*bp_motif.bed')
  //tuple val(sampleId),path ('*.bedpe'),path ("createMotifGCfile_mqc_versions.yml")

  script:
  """
  #Generate BEDPE files
  bedtools bamtobed -bedpe -i $sampleBam | \\
  awk 'OFS = "\t" {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$6-\$2}' | awk '\$11 >=0' > $strBedPE

  #Get GC content
  awk 'OFS = "\t" {print \$1, \$2, \$6, \$7, \$11}' $strBedPE | \\
  sort -k1,1 -k2,2n | \\
  bedtools nuc -fi $genomeFile -bed - | \\
  awk 'OFS = "\t" {print \$1, \$2, \$3, \$4, \$5, \$7}' > $strBed

  #Filter bedfiles for GC content
  awk 'BEGIN {FS=OFS="\t"} FNR==NR {arr[\$4]=\$6;next} (\$7 in arr) {print \$0, arr[\$7]}' \\
  $strBed $strBedPE> $strBedFilterPE

  #Generate 4 mer by collapsing paired reads into 5' nmer and then getting nucleotides from fasta files
  awk -v nmer="${params.nmer}" 'OFS="\t" {print \$1, \$2, \$2+${params.nmer}, \$7, \$8, \$9, \$11, \$12, \$1, \$2, \$6}' $strBedFilterPE > $strBPr1
  awk -v nmer="${params.nmer}" 'OFS="\t" {print \$4, \$6-${params.nmer}, \$6, \$7, \$8, \$10, \$11, \$12, \$1, \$2, \$6}' $strBedFilterPE > $strBPr2
  bedtools getfasta -fi $genomeFile -bed $strBPr1 -s -bedOut -fo | awk 'OFS="\t" {print \$1, \$10, \$11, \$6, \$7, \$8, toupper(\$12)}' - > $strBPr1FA
  bedtools getfasta -fi $genomeFile -bed $strBPr2 -s -bedOut -fo | awk 'OFS="\t" {print toupper(\$12)}' - > $strBPr2FA

  #Put it all together
  awk '{getline line < f2; print \$0 "\t" line}' f2="$strBPr2FA" "$strBPr1FA" > "$strBPmotif"
  
  cat <<-END_VERSIONS > createMotifGCfile_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}