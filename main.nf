 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

include {lib_complex} from './modules/lib_complex'
include {fastqc} from './modules/fastqc'
include {align} from './modules/align'
include {sort_bam} from './modules/sort_bam'
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
include {uropa} from './modules/uropa'
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

    chPeakFiles = peak_bed_graph(chDACFilteredFiles) 
    
    uropa(chPeakFiles,chGeneAnotation) //Verify if it is necessary if its helpful
    chBedFiles = bam_to_bed(chDACFilteredFiles) //
    chUniqueFrags = unique_frags(chBedFiles).collect() //
    chPeakAllFiles = chPeakFiles.collect() //

    
    //Fragments and peaks Plot *******************************************************
    chNarrowPeakFiles = chPeakAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.narrowPeak') }} // Filter the narrowPeak files
    chFragAndPeaksFilesReport = frags_and_peaks(chNarrowPeakFiles,chUniqueFrags,\
    chMultiQCFragsHeader,chMultiQCPeaksHeader,chReportFragPeaks,chSampleInfo)
    //*********************************************************************************

    // SNP Fingerprint and plot process ***************************************************
    /*chSnpFingerprintComplete = snp_fingerprint(chIndexFiles, chSNPS_ref, chGenome).collect() 
    chSnpFingerprintCompleteAllfiles = chSnpFingerprintComplete.collect()
    chVCFGZFiles = chSnpFingerprintCompleteAllfiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.vcf.gz') }} // Filter the vcf files
    chFootPrintPDF = snp_footprint_clustering(chVCFGZFiles,chRSNPFootprint,chSampleInfo)*/

    //ENRICHMENT      ***************************************************
    chEnrichmentFilesCSV = enrichment(chDACFilteredFiles,chEnrichmentScript).collect()
    chEnrichmentFilesReport = enrichmentReport(chSampleInfo,chEnrichmentFilesCSV,chReportEnrichment).collect()
    chMergedEnrichmentReport = merge_enrichment_reports(chEnrichmentFilesReport,chMultiQCEnrichmentHeader,chMergeReportEnrichment,chSampleInfo).collect()
    

    //Final Report

    // Definir o caminho de entrada
    chAllPreviousFiles = Channel.fromPath("${workflow.projectDir}/${params.outputFolder}/")

    chFinalReport = multiqc(chIGVReportMerged,chFragmentsSizeFiles,
        chSNPSMaSHPlot,chEnrichmentFilesReport,chFragAndPeaksFilesReport,chMultiQCConfig,chAllPreviousFiles)

    moveSoftFiles(chFinalReport)
    
}

