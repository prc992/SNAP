 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

include {lib_complex} from './modules/lib_complex'
include {fastqc} from './modules/fastqc'
include {align} from './modules/align'
include {sort_bam} from './modules/sort_bam'
include {sort_readname_bam} from './modules/sort_bam'
include {unique_sam} from './modules/unique_sam'
include {enrichment} from './modules/enrichment'
include {index_sam} from './modules/index_sam'
include {dedup} from './modules/dedup'
include {dac_exclusion} from './modules/dac_exclusion'
include {fetch_chrom_sizes} from './modules/fetch_chrom_sizes'
include {peak_bed_graph} from './modules/peak_bed_graph'
include {bam_to_bed} from './modules/bam_to_bed'
include {unique_frags} from './modules/unique_frags'
include {trim} from './modules/trim'
include {snp_fingerprint} from './modules/snp_fingerprint'
include {bedgraph_to_bigwig} from './modules/bedgraph_to_bigwig'
include {lenght_fragment_dist_step1} from './modules/lenght_fragment_dist_step'
include {lenght_fragment_dist_step2} from './modules/lenght_fragment_dist_step'
include {pileups_report} from './modules/pileups_report'
include {multiqc} from './modules/multiqc'
include {snp_footprint_clustering} from './modules/snp_footprint_clustering'

include {downloadSNPRef} from './modules/download'
include {downloadDACFile} from './modules/download'
include {downloadGeneAnotation} from './modules/download'
include {downloadGenome} from './modules/download'

include {createGenomeIndex} from './modules/createGenomeIndex'
include {createSamplesheet} from './modules/createSamplesheet'

include {createStatsSamtoolsfiltered} from './modules/createStatsSamtoolsfiltered'
include {quality_filter} from './modules/quality_filter'
include {lib_complex_preseq} from './modules/lib_complex_preseq'
include {calcFragsLength} from './modules/calcFragsLength'
include {fragLenHist} from './modules/fragLenHist'
include {frags_and_peaks} from './modules/frags_and_peaks'
include {enrichmentReport} from './modules/enrichmentReport'
include {merge_enrichment_reports} from './modules/merge_enrichment_reports'
include {bam_to_bedgraph} from './modules/bam_to_bedgraph'
include {igv_reports} from './modules/igv_reports'
include {igv_sample_reports} from './modules/igv_reports'
include {igv_consolidate_report} from './modules/igv_reports'
include {igv_session} from './modules/igv_reports'
include {moveSoftFiles} from './modules/moveSoftFiles'

process createSMaSHFingerPrint{
    label 'med_cpu_med_mem'

    container = params.containers.python

    tag "All Samples"   

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/SMaSH/", mode : 'copy'
    
    input:
    path (chSNPSMaSH)
    path (chSNPS_ref)
    path (chBamAndBai)

    output:
    path('pval_out.txt')

    script:
    """
    python3 $chSNPSMaSH -i $chSNPS_ref ALL
    """
}

process createSMaSHFingerPrintPlot{
    label 'med_cpu_med_mem'

    container = params.containers.snp_smash_plot

    tag "All Samples"   

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/SMaSH/", mode : 'copy'
    
    input:
    path (chSMaSHOutout)
    path (chSNPSMaSHPyPlot)

    output:
    path('*.jpg')

    script:
    """
    python3 $chSNPSMaSHPyPlot -i $chSMaSHOutout -o Dendrogram_of_Samples_by_SNP_Profile.jpg
    """
}

process createMotifGCfile {
  label 'process_medium'

  //Docker Image
  container = params.containers.bedtools

  tag "Sample - $sampleId"  

  publishDir "${workflow.projectDir}/${params.outputFolder}/motifs/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId),path(sampleBam),val(_)
  each path (genomeFile)
  each path (genomeIndexFiles)

  exec:
  String strBed = sampleId + '_frags_gc.bed.bed'
  String strBedPE = sampleId + '.bedpe'
  String strBedFilterPE = sampleId + '_filtered.bedpe'
  String nmer = 4 
  String strBPr1 = sampleId + '_' + nmer + '_bp_r1.bed'
  String strBPr2 = sampleId + '_' + nmer + '_bp_r2.bed'
  String strBPr1FA = sampleId + '_' + nmer + 'NMER_bp_r1.fa.bed'
  String strBPr2FA = sampleId + '_' + nmer + 'NMER_bp_r2.fa.bed'
  String strBPmotif = sampleId + '_' + nmer + 'NMER_bp_motif.bed'
  
  output:
  path ('*.*')
  path (sampleBam)
  //tuple val(sampleId),path ('*.bedpe'),path ("createMotifGCfile_mqc_versions.yml")

  script:
  """
  #Generate BEDPE files
  bedtools bamtobed -bedpe -i $sampleBam | \\
  awk 'OFS = "\t" {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$6-\$2}' | awk '\$11 >=0' > $strBedPE

  #Get GC content
  awk 'OFS = "\t" {print \$1, \$2, \$6, \$7, \$11}' $strBedPE | \\
  sort -k1,1 -k2,2n | bedtools nuc -fi $genomeFile -bed - | \\
  awk 'OFS = "\t" {print \$1, \$2, \$3, \$4, \$5, \$7}' > $strBed

  #Filter bedfiles for GC content
  awk 'BEGIN {FS=OFS="\t"} FNR==NR {arr[\$4]=\$6;next} (\$7 in arr) {print \$0, arr[\$7]}' \\
  $strBed $strBedPE> $strBedFilterPE

  #Generate 4 mer by collapsing paired reads into 5' nmer and then getting nucleotides from fasta files
  awk -v nmer="$nmer" 'OFS="\t" {print \$1, \$2, \$2+$nmer, \$7, \$8, \$9, \$11, \$12, \$1, \$2, \$6}' $strBedFilterPE > $strBPr1
  awk -v nmer="$nmer" 'OFS="\t" {print \$4, \$6-$nmer, \$6, \$7, \$8, \$10, \$11, \$12, \$1, \$2, \$6}' $strBedFilterPE > $strBPr2
  bedtools getfasta -fi $genomeFile -bed $strBPr1 -s -bedOut -fo | awk 'OFS="\t" {print \$1, \$10, \$11, \$6, \$7, \$8, toupper(\$12)}' - > $strBPr1FA
  bedtools getfasta -fi $genomeFile -bed $strBPr2 -s -bedOut -fo | awk 'OFS="\t" {print toupper(\$12)}' - > $strBPr2FA

  #Put it all together
  awk '{getline line < f2; print \$0 "\t" line}' f2="$strBPr2FA" "$strBPr1FA" > "$strBPmotif"
  
  cat <<-END_VERSIONS > createMotifGCfile_mqc_versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
  END_VERSIONS
  """

  /*
  #Get GC content
  awk 'OFS = "\t" {print \$1, \$2, \$6, \$7, \$11}' $strBedPE | \\
  sort -k1,1 -k2,2n | bedtools nuc -fi $genomeFile -bed - | \\
  awk 'OFS = "\t" {print \$1, \$2, \$3, \$4, \$5, \$7}' > $strBed

  #Filter bedfiles for GC content
  awk 'BEGIN {FS=OFS="\t"} FNR==NR {arr[\$4]=\$6;next} (\$7 in arr) {print \$0, arr[\$7]}' \\
  $strBed $strBedPE> $strBedFilterPE

  #Generate 4 mer by collapsing paired reads into 5' nmer and then getting nucleotides from fasta files
  awk -v nmer="$nmer" 'OFS="\t" {print \$1, \$2, \$2+$nmer, \$7, \$8, \$9, \$11, \$12, \$1, \$2, \$6}' $strBedFilterPE > $strBPr1
  awk -v nmer="$nmer" 'OFS="\t" {print \$4, \$6-$nmer, \$6, \$7, \$8, \$10, \$11, \$12, \$1, \$2, \$6}' $strBedFilterPE > $strBPr2
  bedtools getfasta -fi $genomeFile -bed $strBPr1 -s -bedOut -fo | awk 'OFS="\t" {print \$1, \$10, \$11, \$6, \$7, \$8, toupper(\$12)}' - > $strBPr1FA
  bedtools getfasta -fi $genomeFile -bed $strBPr2 -s -bedOut -fo | awk 'OFS="\t" {print toupper(\$12)}' - > $strBPr2FA

  #Put it all together
  awk '{getline line < f2; print \$0 "\t" line}' f2="$strBPr2FA" "$strBPr1FA" > "$strBPmotif"
  */
}


workflow {
    // Static information about the pipeline
    def githubPath = "https://github.com/prc992/SNAP"
    def releaseVersion = "v1.0.16"

    // ASCII art for SNAP
    def asciiArt = """
      ███████╗ ███╗   ██╗ █████╗ ██████╗ 
      ██╔════╝ ████╗  ██║██╔══██╗██╔══██╗
      ███████╗ ██╔██╗ ██║███████║██████╔╝
      ╚════██║ ██║╚██╗██║██╔══██║██╔═══╝ 
      ███████╗ ██║ ╚████║██║  ██║██║     
      ╚══════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     
    """

    // Print the introductory message
    println asciiArt
    println "SNAP pipeline running, created by BacaLab. https://bacalab.dana-farber.org/"
    println "SNAP: Streamlined Nextflow Analysis Pipeline for profiling circulating histone modifications identifies tumor epigenomic signatures in cancer plasma."
    println "GitHub repository: ${githubPath}"
    println "Release version: ${releaseVersion}"


    //Auxiliar code
    chEnrichmentScript= Channel.fromPath("$params.pathEnrichmentScript")
    chRfrag_plotFragDist = Channel.fromPath("$params.pathRfrag_plotFragDist")
    chRComparison = Channel.fromPath("$params.pathRComparison")
    chRPileups= Channel.fromPath("$params.pathRPileups")
    //chRSNPFootprint = Channel.fromPath("$params.pathSNPFootprint")
    chSNPSMaSH = Channel.fromPath("$params.pathSNPSMaSH")
    chSNPSMaSHPyPlot = Channel.fromPath("$params.pathSNPSMaSHPlot")
    chReportFragHist = Channel.fromPath("$params.pathReportFragHist")
    chReportFragPeaks = Channel.fromPath("$params.pathReportFragPeaks")
    chReportEnrichment = Channel.fromPath("$params.pathReportEnrichment")
    chMergeReportEnrichment = Channel.fromPath("$params.pathMergeReportEnrichment")
    chIGVFilestoSessions = Channel.fromPath("$params.pathIGVFilestoSessions")

    //Assets
    chPileUpBED = Channel.fromPath("$params.genes_pileup_report")
    chMultiQCConfig = Channel.fromPath("$params.multiqc_config")
    chMultiQCHousekeepingHeader = Channel.fromPath("$params.multiqc_housekeeping_header")
    chMultiQCFragLenHeader = Channel.fromPath("$params.multiqc_frag_len_header")
    chMultiQCPeaksHeader = Channel.fromPath("$params.multiqc_tot_frag_header")
    chMultiQCFragsHeader = Channel.fromPath("$params.multiqc_tot_peaks_header")
    chMultiQCEnrichmentHeader = Channel.fromPath("$params.multiqc_enrichment_header")
    
    // Create the genome directory if it doesn't exist
    """
    mkdir -p ${projectDir}/ref_files/genome
    """.execute().waitFor()
    refDir = Channel.fromPath("${projectDir}/ref_files/genome")

    // Read the GenomePaths spreadsheet and filter the row that matches the genome parameter
    chGenomesSheet = Channel.fromPath(params.genomeInfoPaths)
    chGenomesInfo = chGenomesSheet \
        | splitCsv(header:true) \
        | filter { row -> row.Genome == params.genome } \
        | ifEmpty { error "No matching Genome found in the GenomePaths spreadsheet. Exiting workflow." }
        | map { row-> tuple(row.Genome,row.faGZFile,row.GeneAnotation, row.DACList,row.SNP) }
    // Destructure and store each column into separate variables
    chGenomesInfo
        .map { genome, faGZFile, geneAnnotation, dacList, snp ->
            [genome, faGZFile, geneAnnotation, dacList, snp]
        }

    // Download the genome, gene annotation, and DAC file
    chGenome = downloadGenome(chGenomesInfo,refDir)
    chGenomeIndex = createGenomeIndex(chGenomesInfo,chGenome,refDir)
    chGeneAnotation = downloadGeneAnotation(chGenomesInfo,refDir)
    chChromSizes = fetch_chrom_sizes(chGenomesInfo,refDir)
    chDACFileRef = downloadDACFile(chGenomesInfo,refDir)
    
    // If the 'samplesheet' parameter is provided, use it directly; otherwise, create a new samplesheet
    if (params.samplesheet) {
        //println "Using provided samplesheet: ${params.samplesheet}"
        chSampleSheet = Channel.fromPath(params.samplesheet)
    } else {
        //println "Creating samplesheet because none was provided."
        chSampleSheet = createSamplesheet(
            params.sample_dir, 
            params.enrichment_mark ?: 'no_enrichment_mark'
        )
    }

    // Read the SampleSheet provided by the user or created by the pipeline
    chSampleInfo = chSampleSheet \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.enrichment_mark, row.read1, row.read2) }

    //chSNPS_ref = downloadSNPRef(chGenomesInfo,chSampleInfo)
    chSNPS_ref = downloadSNPRef(chGenomesInfo)

    fastqc(chSampleInfo) 
    chTrimFiles = trim(chSampleInfo)
    chAlignFiles = align(chTrimFiles,chGenome,chGenomeIndex) 
    chSortedFiles = sort_bam(chAlignFiles)
    lib_complex(chSortedFiles) 
    lib_complex_preseq(chSortedFiles) 
    chUniqueFiles = unique_sam(chSortedFiles) 

    chFilteredFiles = quality_filter(chUniqueFiles) 
    chStatsSamtools = createStatsSamtoolsfiltered(chFilteredFiles) 
    chDedupFiles = dedup(chFilteredFiles) 
    chDACFilteredFiles = dac_exclusion(chDedupFiles,chDACFileRef) 
    chIndexFiles = index_sam(chDACFilteredFiles)

    
    chAllIndexFiles = chIndexFiles.collect()
    chAllBAMandBAIIndexFiles = chAllIndexFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }}

    chSMaSHOutout = createSMaSHFingerPrint(chSNPSMaSH,chSNPS_ref,chAllBAMandBAIIndexFiles)
    chSNPSMaSHPlot = createSMaSHFingerPrintPlot(chSMaSHOutout,chSNPSMaSHPyPlot)

    chBedGraphFiles = bam_to_bedgraph(chIndexFiles)
    chBigWig = bedgraph_to_bigwig(chBedGraphFiles,chChromSizes)

    
    //Pileups ****************************************************************
    chIGVReportsHtml = igv_sample_reports(chBedGraphFiles,chPileUpBED,chGenome,chGenomeIndex).collect()
    chIGVReportMerged = igv_consolidate_report(chIGVReportsHtml,chMultiQCHousekeepingHeader)

    chBigWigAllFiles = chBigWig.collect()
    chBigWigOnlyFiles = chBigWigAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bw') }} // Filter the tdf files
    chIGVSession = igv_session(chBigWigOnlyFiles,chIGVFilestoSessions,chGenomesInfo,chPileUpBED)
        
    
    //Fragment Length Distribution ************************************************
    chFragmentsSize = calcFragsLength(chIndexFiles).collect()
    chFragmentsSizeFiles = chFragmentsSize.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.fragment_sizes.txt') }} // Filter the Fragments Size files
    //************************************************************************

    // Create a channel with a tuple containing a string ID, the BAM file, and an additional string
    chBamTest = Channel.fromPath('/Users/prc992/Desktop/DFCI/2-SNAP/9-MotifGunTest/HS_cK20_AM_MH_unique_sorted_deduped_filtered.bam')
                    .map { file -> tuple('sample_test', file, 'alo') }

    chNameSortedFiles = sort_readname_bam(chBamTest)
    createMotifGCfile(chNameSortedFiles, chGenome, chGenomeIndex)
    //createMotifGCfile(chDACFilteredFiles, chGenome, chGenomeIndex)

    chPeakFiles = peak_bed_graph(chDACFilteredFiles) 
    chBedFiles = bam_to_bed(chDACFilteredFiles) //
    
    /*chUniqueFrags = unique_frags(chBedFiles).collect() //
    chPeakAllFiles = chPeakFiles.collect() //

    
    //Fragments and peaks Plot *******************************************************
    chNarrowPeakFiles = chPeakAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.narrowPeak') }} // Filter the narrowPeak files
    chFragAndPeaksFilesReport = frags_and_peaks(chNarrowPeakFiles,chUniqueFrags,\
    chMultiQCFragsHeader,chMultiQCPeaksHeader,chReportFragPeaks,chSampleInfo)
    //*********************************************************************************

    //ENRICHMENT      ***************************************************
    chEnrichmentFilesCSV = enrichment(chDACFilteredFiles,chEnrichmentScript).collect()
    chEnrichmentFilesReport = enrichmentReport(chSampleInfo,chEnrichmentFilesCSV,chReportEnrichment).collect()
    chMergedEnrichmentReport = merge_enrichment_reports(chEnrichmentFilesReport,chMultiQCEnrichmentHeader,chMergeReportEnrichment,chSampleInfo).collect()
    

    //Final Report

    // Definir o caminho de entrada
    chAllPreviousFiles = Channel.fromPath("${workflow.projectDir}/${params.outputFolder}/")

    chFinalReport = multiqc(chIGVReportMerged,chFragmentsSizeFiles,
        chSNPSMaSHPlot,chEnrichmentFilesReport,chFragAndPeaksFilesReport,chMultiQCConfig,chAllPreviousFiles)

    moveSoftFiles(chFinalReport)*/
    
}

