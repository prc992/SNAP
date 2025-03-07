 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

//local modules
include {moveSoftFiles} from './modules/moveSoftFiles'

//subworkflows
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
    chMultiQCFragsHeader = Channel.fromPath("$params.multiqc_tot_frag_header")
    chMultiQCPeaksHeader = Channel.fromPath("$params.multiqc_tot_peaks_header")
    chMultiQCEnrichmentHeader = Channel.fromPath("$params.multiqc_enrichment_header")

    def steps = ['INITIALIZATION', 'DOWNLOAD_REFERENCES', 'BAM_PROCESSING', 'BAM_SIGNAL_PROCESSING', 'FRAGMENTS_PROCESSING']
    def run_steps = steps.takeWhile { it != params.until } + params.until


    if ('INITIALIZATION' in run_steps) {
        INITIALIZATION()
        chGenomesInfo = INITIALIZATION.out.genomes_info
        refDir = INITIALIZATION.out.ref_dir
        chSampleInfo = INITIALIZATION.out.sample_info
        chFastaQC = INITIALIZATION.out.fastqc_files
        chFilesReportInitialization = INITIALIZATION.out.files_report_initialization
        chInitReport = INITIALIZATION.out.init_report
        }

    if ('DOWNLOAD_REFERENCES' in run_steps) {
        DOWNLOAD_REFERENCES(chGenomesInfo,refDir)

        chGenome = DOWNLOAD_REFERENCES.out.genome
        chGenomeIndex = DOWNLOAD_REFERENCES.out.genome_index
        chChromSizes = DOWNLOAD_REFERENCES.out.chrom_sizes
        chDACFileRef = DOWNLOAD_REFERENCES.out.dac_file_ref
        chSampleInfo = DOWNLOAD_REFERENCES.out.sample_info
        chSNPS_ref = DOWNLOAD_REFERENCES.out.snp_ref
        }

    if ('BAM_PROCESSING' in run_steps) {
        BAM_PROCESSING (chSampleInfo, chGenome, chGenomeIndex,chChromSizes,chDACFileRef,chSNPSMaSH,chSNPS_ref,chSNPSMaSHPyPlot,\
                        chFilesReportInitialization,chMultiQCConfig,chInitReport)

        chBAMProcessedFiles = BAM_PROCESSING.out.bam_processed
        chBAMProcessedIndexFiles = BAM_PROCESSING.out.bam_processed_index
        chSNPSMaSHPlot = BAM_PROCESSING.out.report_SNP_SMaSH
        chLibComplexPreseq = BAM_PROCESSING.out.lib_complex
        chFilesReportBamProcessing = BAM_PROCESSING.out.files_report_bam_processing
        chBAMProcessReport = BAM_PROCESSING.out.bam_process_report
        }

    if ('BAM_SIGNAL_PROCESSING' in run_steps) {
        BAM_SIGNAL_PROCESSING(chSampleInfo,chBAMProcessedFiles,chBAMProcessedIndexFiles,chChromSizes,chPileUpBED,chGenome,chGenomeIndex,\
                            chMultiQCHousekeepingHeader,chIGVFilestoSessions,chGenomesInfo,chMultiQCPeaksHeader,chReportPeaks,\
                            chEnrichmentScript,chReportEnrichment,chMergeReportEnrichment,chMultiQCEnrichmentHeader,\
                            chFilesReportBamProcessing,chFilesReportInitialization,chMultiQCConfig,chBAMProcessReport)

        chIGVReportMerged = BAM_SIGNAL_PROCESSING.out.igv_report_merged
        chEnrichmentFilesReport = BAM_SIGNAL_PROCESSING.out.merge_enrichment_reports
        chPeaksReport = BAM_SIGNAL_PROCESSING.out.peaks_report
        chFilesReportSignalProcess = BAM_SIGNAL_PROCESSING.out.files_report_bam_signal_processing
        chBAMSignalReport = BAM_SIGNAL_PROCESSING.out.bam_signal_report
        }

    if ('FRAGMENTS_PROCESSING' in run_steps) {
        FRAGMENTS_PROCESSING(chBAMProcessedFiles,chBAMProcessedIndexFiles,chGenome,chGenomeIndex,\
                            chMultiQCFragsHeader,chReportFrags,chFilesReportSignalProcess,\
                            chFilesReportBamProcessing,chFilesReportInitialization,chMultiQCConfig,chBAMSignalReport)

        chFragmentsSizeFiles = FRAGMENTS_PROCESSING.out.frag_size_files
        chFragReport = FRAGMENTS_PROCESSING.out.frag_report
        }

    /*// Criamos um canal que só será ativado quando todas as saídas estiverem prontas
    chFinalTrigger = chIGVReportMerged
        .combine(chFastaQC)
        .combine(chFragmentsSizeFiles)
        .combine(chSNPSMaSHPlot)
        .combine(chLibComplexPreseq)
        .combine(chEnrichmentFilesReport)
        .combine(chPeaksReport)
        .combine(chFragReport)

        chAllPreviousFiles = Channel.fromPath("${workflow.projectDir}/${params.outputFolder}/")
        chFinalReport = multiqc(chFinalTrigger, chMultiQCConfig, chAllPreviousFiles)
        moveSoftFiles(chFinalReport)*/
    
}

