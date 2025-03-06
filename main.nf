 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

//local modules
include {multiqc} from './modules/multiqc'
include {moveSoftFiles} from './modules/moveSoftFiles'

//subworkflows
include { INITIALIZATION } from './subworkflows/local/initialization'
include { DOWNLOAD_REFERENCES } from './subworkflows/local/download_references'
include { BAM_PROCESSING } from './subworkflows/local/bam_processing'
include { BAM_SIGNAL_PROCESSING } from './subworkflows/local/bam_signal_process'
include { FRAGMENTS_PROCESSING } from './subworkflows/local/fragments_processing'

process CREATE_DUMMY_FILE {

    input:
    val (file_name)

    output:
    path "*.txt"

    script:
    """
    touch "$file_name".txt
    """
}

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
    chMultiQCDummyFile = Channel.fromPath("$params.multiqc_dummy_file")

    def steps = ['INITIALIZATION', 'DOWNLOAD_REFERENCES', 'BAM_PROCESSING', 'BAM_SIGNAL_PROCESSING', 'FRAGMENTS_PROCESSING']
    def run_steps = steps.takeWhile { it != params.until } + params.until

    // Criamos um canal que serve como gatilho para o processo final
    chFinalTrigger = Channel.empty()

    // Inicializa os canais vazios para evitar erro de variável não definida
    chIGVReportMerged = Channel.empty()
    chFragmentsSizeFiles = Channel.empty()
    chSNPSMaSHPlot = Channel.empty()
    chEnrichmentFilesReport = Channel.empty()
    chPeaksReport = Channel.empty()
    chFragReport = Channel.empty()
    chFastaQC = Channel.empty()
    chLibComplexPreseq = Channel.empty()

    if ('INITIALIZATION' in run_steps) {
        INITIALIZATION()
        chGenomesInfo = INITIALIZATION.out.genomes_info
        refDir = INITIALIZATION.out.ref_dir
        chSampleInfo = INITIALIZATION.out.sample_info
        chFastaQC = INITIALIZATION.out.fastqc_files
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
        BAM_PROCESSING (chSampleInfo, chGenome, chGenomeIndex,chChromSizes,chDACFileRef,chSNPSMaSH,chSNPS_ref,chSNPSMaSHPyPlot)

        chBAMProcessedFiles = BAM_PROCESSING.out.bam_processed
        chBAMProcessedIndexFiles = BAM_PROCESSING.out.bam_processed_index
        chSNPSMaSHPlot = BAM_PROCESSING.out.report_SNP_SMaSH
        chLibComplexPreseq = BAM_PROCESSING.out.lib_complex
        }

    if ('BAM_SIGNAL_PROCESSING' in run_steps) {
        BAM_SIGNAL_PROCESSING(chSampleInfo,chBAMProcessedFiles,chBAMProcessedIndexFiles,chChromSizes,chPileUpBED,chGenome,chGenomeIndex,\
                            chMultiQCHousekeepingHeader,chIGVFilestoSessions,chGenomesInfo,chMultiQCPeaksHeader,chReportPeaks,\
                            chEnrichmentScript,chReportEnrichment,chMergeReportEnrichment,chMultiQCEnrichmentHeader)

        chIGVReportMerged = BAM_SIGNAL_PROCESSING.out.igv_report_merged
        chEnrichmentFilesReport = BAM_SIGNAL_PROCESSING.out.merge_enrichment_reports
        chPeaksReport = BAM_SIGNAL_PROCESSING.out.peaks_report
        
        }

    if ('FRAGMENTS_PROCESSING' in run_steps) {
        FRAGMENTS_PROCESSING(chBAMProcessedFiles,chBAMProcessedIndexFiles,chGenome,chGenomeIndex,\
                            chMultiQCFragsHeader,chReportFrags)

        chFragmentsSizeFiles = FRAGMENTS_PROCESSING.out.frag_size_files
        chFragReport = FRAGMENTS_PROCESSING.out.frag_report
        }

    //
    // Criamos arquivos dummy para evitar bloqueio do `combine()`
    chIGVReportMerged = chIGVReportMerged.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_igv.txt"))
    chFastaQC = chFastaQC.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_fasta.txt"))
    chFragmentsSizeFiles = chFragmentsSizeFiles.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_frag.txt"))
    chSNPSMaSHPlot = chSNPSMaSHPlot.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_snps.txt"))
    chLibComplexPreseq = chLibComplexPreseq.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_libcomplex.txt"))
    chEnrichmentFilesReport = chEnrichmentFilesReport.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_enrichment.txt"))
    chPeaksReport = chPeaksReport.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_peaks.txt"))
    chFragReport = chFragReport.ifEmpty(Channel.fromPath("${workflow.projectDir}/ref_files/multiqc/dummy_frag_report.txt"))

    // Criamos um canal que só será ativado quando todas as saídas estiverem prontas
    // Criamos um canal que só será ativado quando todas as saídas estiverem prontas
    chFinalTrigger = chIGVReportMerged
        .combine(chFastaQC)
        .combine(chFragmentsSizeFiles)
        .combine(chSNPSMaSHPlot)
        .combine(chLibComplexPreseq)
        .combine(chEnrichmentFilesReport)
        .combine(chPeaksReport)
        .combine(chFragReport)
    //

    /*

    //Final Report
    chAllPreviousFiles = Channel.fromPath("${workflow.projectDir}/${params.outputFolder}/")
    chFinalTrigger = chIGVReportMerged
        .mix(chFastaQC)
        .mix(chFragmentsSizeFiles)
        .mix(chSNPSMaSHPlot)
        .mix(chLibComplexPreseq)
        .mix(chEnrichmentFilesReport)
        .mix(chPeaksReport)
        .mix(chFragReport)

    chFinalReport = multiqc(chFinalTrigger,chMultiQCConfig,chAllPreviousFiles)


    /*chFinalReport = multiqc(chIGVReportMerged,chFragmentsSizeFiles,
        chSNPSMaSHPlot,chEnrichmentFilesReport,chPeaksReport,chFragReport,chMultiQCConfig,chAllPreviousFiles)*/

    /*moveSoftFiles(chFinalReport)*/
    
}

