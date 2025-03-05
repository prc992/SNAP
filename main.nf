 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

//include {fastqc} from './modules/fastqc'


//include {sort_readname_bam} from './modules/sort_bam'

//include {enrichment} from './modules/enrichment'


//include {bam_to_bed} from './modules/bam_to_bed'
//include {unique_frags} from './modules/unique_frags'

/// BAM PROCESSING
//include {trim} from './modules/trim'
//include {align} from './modules/align'
//include {sort_bam} from './modules/sort_bam'
//include {lib_complex_preseq} from './modules/lib_complex_preseq'
//include {unique_sam} from './modules/unique_sam'
//include {quality_filter} from './modules/quality_filter'
//include {createStatsSamtoolsfiltered} from './modules/createStatsSamtoolsfiltered'
//include {dedup} from './modules/dedup'
//include {dac_exclusion} from './modules/dac_exclusion'


include {lenght_fragment_dist_step1} from './modules/lenght_fragment_dist_step'
include {lenght_fragment_dist_step2} from './modules/lenght_fragment_dist_step'
include {pileups_report} from './modules/pileups_report'
include {multiqc} from './modules/multiqc'

/// DOWNLOAD_REFERENCES
//include {downloadSNPRef} from './modules/download'
//include {downloadDACFile} from './modules/download'
//include {downloadGeneAnotation} from './modules/download'
//include {downloadGenome} from './modules/download'
//include {createGenomeIndex} from './modules/createGenomeIndex'
//include {createSamplesheet} from './modules/createSamplesheet'
//include {fetch_chrom_sizes} from './modules/fetch_chrom_sizes'

/// BAM SIGNAL PROCESSING
//include {index_sam} from './modules/index_sam'
//include {bam_to_bedgraph} from './modules/bam_to_bedgraph'
//include {bedgraph_to_bigwig} from './modules/bedgraph_to_bigwig'
//include {igv_reports} from './modules/igv_reports'
//include {igv_sample_reports} from './modules/igv_reports'
//include {igv_consolidate_report} from './modules/igv_reports'
//include {igv_session} from './modules/igv_reports'
//include {peak_bed_graph} from './modules/peak_bed_graph'


//include {calcFragsLengthDistribuition} from './modules/calcFragsLength'
include {fragLenHist} from './modules/fragLenHist'
//include {frags_and_peaks} from './modules/frags_and_peaks' //remove definitely
//include {enrichmentReport} from './modules/enrichmentReport'
//include {merge_enrichment_reports} from './modules/merge_enrichment_reports'


include {moveSoftFiles} from './modules/moveSoftFiles'
//include {createSMaSHFingerPrint} from './modules/snp_smash_fingerprint'
//include {createSMaSHFingerPrintPlot} from './modules/snp_smash_fingerprint'
//include {createMotifGCfile} from './modules/end_motif_gc'

include { INITIALIZATION } from './subworkflows/local/initialization'
include { DOWNLOAD_REFERENCES } from './subworkflows/local/download_references'
include { BAM_PROCESSING } from './subworkflows/local/bam_processing'
include { BAM_SIGNAL_PROCESSING } from './subworkflows/local/bam_signal_process'
include { FRAGMENTS_PROCESSING } from './subworkflows/local/fragments_processing'


workflow {
    // Static information about the pipeline
    def githubPath = "https://github.com/prc992/SNAP"
    def releaseVersion = "v1.2.0 - All processes in the pipeline"

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
    chSNPSMaSH = Channel.fromPath("$params.pathSNPSMaSH")
    chSNPSMaSHPyPlot = Channel.fromPath("$params.pathSNPSMaSHPlot")
    chReportFragHist = Channel.fromPath("$params.pathReportFragHist")
    chReportFrags = Channel.fromPath("$params.pathReportFrags")
    chReportPeaks = Channel.fromPath("$params.pathReportPeaks")
    chReportEnrichment = Channel.fromPath("$params.pathReportEnrichment")
    chMergeReportEnrichment = Channel.fromPath("$params.pathMergeReportEnrichment")
    chIGVFilestoSessions = Channel.fromPath("$params.pathIGVFilestoSessions")

    //Assets
    chPileUpBED = Channel.fromPath("$params.genes_pileup_report")
    chMultiQCConfig = Channel.fromPath("$params.multiqc_config")
    chMultiQCHousekeepingHeader = Channel.fromPath("$params.multiqc_housekeeping_header")
    //chMultiQCFragLenHeader = Channel.fromPath("$params.multiqc_frag_len_header") remove definitely
    chMultiQCFragsHeader = Channel.fromPath("$params.multiqc_tot_frag_header")
    chMultiQCPeaksHeader = Channel.fromPath("$params.multiqc_tot_peaks_header")
    chMultiQCEnrichmentHeader = Channel.fromPath("$params.multiqc_enrichment_header")

    // Create the samplesheet, run FastQC, gather the genome information
    INITIALIZATION()

    chGenomesInfo = INITIALIZATION.out.genomes_info
    refDir = INITIALIZATION.out.ref_dir
    chSampleInfo = INITIALIZATION.out.sample_info
   
    // Download the genome, gene annotation, and DAC file
    DOWNLOAD_REFERENCES(chGenomesInfo,refDir)

    chGenome = DOWNLOAD_REFERENCES.out.genome
    chGenomeIndex = DOWNLOAD_REFERENCES.out.genome_index
    chChromSizes = DOWNLOAD_REFERENCES.out.chrom_sizes
    chDACFileRef = DOWNLOAD_REFERENCES.out.dac_file_ref
    chSampleInfo = DOWNLOAD_REFERENCES.out.sample_info
    chSNPS_ref = DOWNLOAD_REFERENCES.out.snp_ref


    // Process the BAM files
    BAM_PROCESSING (chSampleInfo, chGenome, chGenomeIndex,chChromSizes,chDACFileRef,chSNPSMaSH,chSNPS_ref,chSNPSMaSHPyPlot)

    chBAMProcessedFiles = BAM_PROCESSING.out.bam_processed
    chBAMProcessedIndexFiles = BAM_PROCESSING.out.bam_processed_index
    chSNPSMaSHPlot = BAM_PROCESSING.out.report_SNP_SMaSH

    // Process the BAM signal
    BAM_SIGNAL_PROCESSING(chSampleInfo,chBAMProcessedFiles,chBAMProcessedIndexFiles,chChromSizes,chPileUpBED,chGenome,chGenomeIndex,\
                            chMultiQCHousekeepingHeader,chIGVFilestoSessions,chGenomesInfo,chMultiQCPeaksHeader,chReportPeaks,\
                            chEnrichmentScript,chReportEnrichment,chMergeReportEnrichment,chMultiQCEnrichmentHeader)

    chIGVReportMerged = BAM_SIGNAL_PROCESSING.out.igv_report_merged
    chEnrichmentFilesReport = BAM_SIGNAL_PROCESSING.out.merge_enrichment_reports
    chPeaksReport = BAM_SIGNAL_PROCESSING.out.peaks_report


    // Process the fragments
    FRAGMENTS_PROCESSING(chBAMProcessedFiles,chBAMProcessedIndexFiles,chGenome,chGenomeIndex,\
                            chMultiQCFragsHeader,chReportFrags)

    chFragmentsSizeFiles = FRAGMENTS_PROCESSING.out.frag_size_files
    chFragReport = FRAGMENTS_PROCESSING.out.frag_report


    //Final Report
    chAllPreviousFiles = Channel.fromPath("${workflow.projectDir}/${params.outputFolder}/")

    chFinalReport = multiqc(chIGVReportMerged,chFragmentsSizeFiles,
        chSNPSMaSHPlot,chEnrichmentFilesReport,chPeaksReport,chFragReport,chMultiQCConfig,chAllPreviousFiles)

    moveSoftFiles(chFinalReport)
    
}

