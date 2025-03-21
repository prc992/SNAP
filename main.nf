 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

//local modules


//subworkflows
include { PREPROCESSING } from './subworkflows/local/preprocessing.nf'
include { DOWNLOAD_REFERENCES } from './subworkflows/local/download_references'
include { ALIGNMENT } from './subworkflows/local/alignment'
include { BAM_PROCESSING } from './subworkflows/local/bam_processing'
include { FRAGMENTS_PROCESSING } from './subworkflows/local/fragments_processing'
include { BAM_SIGNAL_PROCESSING } from './subworkflows/local/bam_signal_process'


workflow  {
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
    println "Running pipeline with profile: ${workflow.profile}"

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

    if (params.samplesheetBams || params.sample_dir_bam) {
        skip_alignment = true
    } else {
        skip_alignment = false
    }

    def steps = ['PREPROCESSING', 'DOWNLOAD_REFERENCES','ALIGNMENT', 'BAM_PROCESSING', 'FRAGMENTS_PROCESSING','BAM_SIGNAL_PROCESSING']
    def run_steps = steps.takeWhile { it != params.until } + params.until
    
    if ('PREPROCESSING' in run_steps) {
        PREPROCESSING(chMultiQCConfig,skip_alignment)
        chGenomesInfo = PREPROCESSING.out.genomes_info
        refDir = PREPROCESSING.out.ref_dir
        chSampleInfo = PREPROCESSING.out.sample_info
        chFastaQC = PREPROCESSING.out.fastqc_files
        chFilesReportInitialization = PREPROCESSING.out.files_report_initialization
        chInitReport = PREPROCESSING.out.init_report
        }

    if ('DOWNLOAD_REFERENCES' in run_steps) {
        DOWNLOAD_REFERENCES(chGenomesInfo,refDir)

        chGenome = DOWNLOAD_REFERENCES.out.genome
        chGenomeIndex = DOWNLOAD_REFERENCES.out.genome_index
        chChromSizes = DOWNLOAD_REFERENCES.out.chrom_sizes
        chDACFileRef = DOWNLOAD_REFERENCES.out.dac_file_ref
        chSNPS_ref = DOWNLOAD_REFERENCES.out.snp_ref
        }

    
    if ('ALIGNMENT' in run_steps) {

        if (skip_alignment) {
            chAlign = chSampleInfo.map { sampleId, enrichment_mark, bam -> 
                        tuple(sampleId, file(bam), enrichment_mark) 
                        }

            chFilesReportAlignment = Channel.of("NO_DATA")
            chAlignmentReport = Channel.of("NO_DATA")

        } else {
            ALIGNMENT (chSampleInfo,chGenome,chGenomeIndex,\
                    chFilesReportInitialization,chInitReport,\
                    chMultiQCConfig)
            chAlign = ALIGNMENT.out.align
            chFilesReportAlignment = ALIGNMENT.out.files_report_alignment
            chAlignmentReport = ALIGNMENT.out.aligment_report
        }
        

        
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
    
    if ('FRAGMENTS_PROCESSING' in run_steps) {

        FRAGMENTS_PROCESSING(chGenome,chGenomeIndex,
                            chBAMProcessedFiles,chBAMProcessedIndexFiles,
                            chMultiQCFragsHeader,chReportFrags,
                            chMultiQCConfig,
                            chFilesReportInitialization,chFilesReportBamProcessing,
                            chBAMProcessReport)

        chFragmentsSizeFiles = FRAGMENTS_PROCESSING.out.frag_size_files
        chFragReport = FRAGMENTS_PROCESSING.out.frag_report
        chFragsProcessReport = FRAGMENTS_PROCESSING.out.frag_process_report
        chFilesReportFragmentsProcess = FRAGMENTS_PROCESSING.out.files_report_fragments_processing
        }

    if ('BAM_SIGNAL_PROCESSING' in run_steps) {
        BAM_SIGNAL_PROCESSING (chSampleInfo, chGenome, chGenomeIndex,chChromSizes,skip_alignment,
                            chBAMProcessedFiles,chBAMProcessedIndexFiles,
                            chGenomesInfo,chMultiQCHousekeepingHeader,chMultiQCEnrichmentHeader,chMultiQCPeaksHeader,chIGVFilestoSessions,
                            chMultiQCConfig,chEnrichmentScript,chPileUpBED,chReportPeaks,chReportEnrichment,chMergeReportEnrichment,
                            chFilesReportInitialization,chFilesReportBamProcessing,chFilesReportFragmentsProcess,
                            chFragReport)



        chIGVReportMerged = BAM_SIGNAL_PROCESSING.out.igv_report_merged
        chEnrichmentFilesReport = BAM_SIGNAL_PROCESSING.out.merge_enrichment_reports
        chPeaksReport = BAM_SIGNAL_PROCESSING.out.peaks_report
        chFilesReportSignalProcess = BAM_SIGNAL_PROCESSING.out.files_report_bam_signal_processing
        chBAMSignalReport = BAM_SIGNAL_PROCESSING.out.bam_signal_report
        }
    
}

