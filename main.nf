 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

//local modules


//subworkflows
include { INITIALIZATION } from './subworkflows/local/initialization'
include { DOWNLOAD_REFERENCES } from './subworkflows/local/download_references'
include { ALIGNMENT } from './subworkflows/local/alignment'
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
    chMultiQCFragsHeader = Channel.fromPath("$params.multiqc_tot_frag_header")
    chMultiQCPeaksHeader = Channel.fromPath("$params.multiqc_tot_peaks_header")
    chMultiQCEnrichmentHeader = Channel.fromPath("$params.multiqc_enrichment_header")

    def steps = ['INITIALIZATION', 'DOWNLOAD_REFERENCES','ALIGNMENT', 'BAM_PROCESSING', 'BAM_SIGNAL_PROCESSING', 'FRAGMENTS_PROCESSING']
    def run_steps = steps.takeWhile { it != params.until } + params.until


    if ('INITIALIZATION' in run_steps) {
        INITIALIZATION(chMultiQCConfig)
        chGenomesInfo = INITIALIZATION.out.genomes_info
        refDir = INITIALIZATION.out.ref_dir
        chSampleInfo = INITIALIZATION.out.sample_info
        chFastaQC = INITIALIZATION.out.fastqc_files
        chFilesReportInitialization = INITIALIZATION.out.files_report_initialization
        chInitReport = INITIALIZATION.out.init_report
        }

    if ('DOWNLOAD_REFERENCES' in run_steps) {
        DOWNLOAD_REFERENCES(chGenomesInfo,refDir)

        chSampleInfo = DOWNLOAD_REFERENCES.out.sample_info
        chGenome = DOWNLOAD_REFERENCES.out.genome
        chGenomeIndex = DOWNLOAD_REFERENCES.out.genome_index
        chChromSizes = DOWNLOAD_REFERENCES.out.chrom_sizes
        chDACFileRef = DOWNLOAD_REFERENCES.out.dac_file_ref
        chSNPS_ref = DOWNLOAD_REFERENCES.out.snp_ref
        }

    
    if ('ALIGNMENT' in run_steps) {
        ALIGNMENT (chSampleInfo,chGenome,chGenomeIndex,\
                    chFilesReportInitialization,chInitReport,\
                    chMultiQCConfig)

        chAlign = ALIGNMENT.out.align
        chFilesReportAlignment = ALIGNMENT.out.files_report_alignment
        chAlignmentReport = ALIGNMENT.out.aligment_report
        }

    if ('BAM_PROCESSING' in run_steps) {
        BAM_PROCESSING (chSampleInfo, chGenome, chGenomeIndex,chChromSizes,chDACFileRef,chSNPS_ref,
                        chAlign,
                        chSNPSMaSH,chSNPSMaSHPyPlot,chMultiQCConfig,
                        chFilesReportInitialization,chInitReport,
                        chFilesReportAlignment,chAlignmentReport)

        chBAMProcessedFiles = BAM_PROCESSING.out.bam_processed
        chBAMProcessedIndexFiles = BAM_PROCESSING.out.bam_processed_index
        chSNPSMaSHPlot = BAM_PROCESSING.out.report_SNP_SMaSH
        chLibComplexPreseq = BAM_PROCESSING.out.lib_complex
        chFilesReportBamProcessing = BAM_PROCESSING.out.files_report_bam_processing
        chBAMProcessReport = BAM_PROCESSING.out.bam_process_report
        }

    if ('BAM_SIGNAL_PROCESSING' in run_steps) {
        BAM_SIGNAL_PROCESSING (chSampleInfo, chGenome, chGenomeIndex,chChromSizes,
                            chBAMProcessedFiles,chBAMProcessedIndexFiles,
                            chGenomesInfo,chMultiQCHousekeepingHeader,chMultiQCEnrichmentHeader,chMultiQCPeaksHeader,chIGVFilestoSessions,
                            chMultiQCConfig,chEnrichmentScript,chPileUpBED,chReportPeaks,chReportEnrichment,chMergeReportEnrichment,
                            chFilesReportInitialization,chFilesReportBamProcessing,chBAMProcessReport)



        chIGVReportMerged = BAM_SIGNAL_PROCESSING.out.igv_report_merged
        chEnrichmentFilesReport = BAM_SIGNAL_PROCESSING.out.merge_enrichment_reports
        chPeaksReport = BAM_SIGNAL_PROCESSING.out.peaks_report
        chFilesReportSignalProcess = BAM_SIGNAL_PROCESSING.out.files_report_bam_signal_processing
        chBAMSignalReport = BAM_SIGNAL_PROCESSING.out.bam_signal_report
        }
/*
    if ('FRAGMENTS_PROCESSING' in run_steps) {
        FRAGMENTS_PROCESSING(chBAMProcessedFiles,chBAMProcessedIndexFiles,chGenome,chGenomeIndex,\
                            chMultiQCFragsHeader,chReportFrags,chFilesReportSignalProcess,\
                            chFilesReportBamProcessing,chFilesReportInitialization,chMultiQCConfig,chBAMSignalReport)

        chFragmentsSizeFiles = FRAGMENTS_PROCESSING.out.frag_size_files
        chFragReport = FRAGMENTS_PROCESSING.out.frag_report
        chFragsProcessReport = FRAGMENTS_PROCESSING.out.frag_process_report
        }

*/
    
}

