process createMotifGCfile {
  label 'process_medium'

  //Docker Image
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/motifs/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),val(enrichment_mark),val(control),val(read_method),path(sampleBam),val (_)
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
  tuple val(sampleId),val(enrichment_mark),val(control),val(read_method),path(sampleBam),path ('*bp_motif.bed'),path ("createMotifGCfile_mqc_versions.yml")

  script:
  def motifCommand = ""
  if (read_method == "PE") {
    motifCommand = """
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
    """
  } else {
    motifCommand = """
    #  Single-End: BED files
    bedtools bamtobed -i $sampleBam > ${sampleId}.bed

    # Expand 5' according to strand
    awk -v nmer="${params.nmer}" 'BEGIN{OFS="\\t"} {
      if (\$6 == "+") {
        print \$1, \$2, \$2 + nmer, \$4, \$5, \$6
      } else {
        print \$1, \$3 - nmer, \$3, \$4, \$5, \$6
      }
    }' ${sampleId}.bed > ${sampleId}_expanded.bed

    # GC content
    awk 'OFS = "\\t" {print \$1, \$2, \$3, \$4, \$5}' ${sampleId}_expanded.bed | \\
    sort -k1,1 -k2,2n | \\
    bedtools nuc -fi $genomeFile -bed - | \\
    awk 'OFS = "\\t" {print \$1, \$2, \$3, \$4, \$5, \$7}' > $strBed
    
    # Extract sequences from the genome
    bedtools getfasta -fi $genomeFile -bed ${sampleId}_expanded.bed -s -bedOut -fo | \\
    awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, toupper(\$8), toupper(\$8)}' > "$strBPmotif"
    """
  }
  
  """
  echo "Generating motif GC file for sample $sampleId ($read_method)"
  $motifCommand 
  
  cat <<-END_VERSIONS > createMotifGCfile_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """
}