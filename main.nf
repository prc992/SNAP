 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

include {fastqc} from './modules/fastqc'


include {sort_readname_bam} from './modules/sort_bam'

include {enrichment} from './modules/enrichment'
include {index_sam} from './modules/index_sam'

include {dac_exclusion} from './modules/dac_exclusion'
include {peak_bed_graph} from './modules/peak_bed_graph'
include {bam_to_bed} from './modules/bam_to_bed'
include {unique_frags} from './modules/unique_frags'

/// BAM PROCESSING
//include {trim} from './modules/trim'
//include {align} from './modules/align'
//include {sort_bam} from './modules/sort_bam'
//include {lib_complex_preseq} from './modules/lib_complex_preseq'
//include {unique_sam} from './modules/unique_sam'
//include {quality_filter} from './modules/quality_filter'
//include {createStatsSamtoolsfiltered} from './modules/createStatsSamtoolsfiltered'
//include {dedup} from './modules/dedup'

include {bedgraph_to_bigwig} from './modules/bedgraph_to_bigwig'
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
include {createSMaSHFingerPrint} from './modules/snp_smash_fingerprint'
include {createSMaSHFingerPrintPlot} from './modules/snp_smash_fingerprint'
include {createMotifGCfile} from './modules/end_motif_gc'

include { DOWNLOAD_REFERENCES } from './subworkflows/local/download_references'
include { BAM_PROCESSING } from './subworkflows/local/bam_processing'

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
    DOWNLOAD_REFERENCES(chGenomesInfo,refDir)

    chGenome = DOWNLOAD_REFERENCES.out.genome
    chGenomeIndex = DOWNLOAD_REFERENCES.out.genome_index
    chChromSizes = DOWNLOAD_REFERENCES.out.chrom_sizes
    chDACFileRef = DOWNLOAD_REFERENCES.out.dac_file_ref
    chSampleInfo = DOWNLOAD_REFERENCES.out.sample_info
    chSNPS_ref = DOWNLOAD_REFERENCES.out.snp_ref

    fastqc(chSampleInfo)

    // Process the BAM files
    BAM_PROCESSING (chSampleInfo, chGenome, chGenomeIndex)

    chDedupFiles = BAM_PROCESSING.out.bam_deduped

    //************************************************************************
    //DOWNLOAD_REFERENCES
    //************************************************************************

    //chGeneAnotation = downloadGeneAnotation(chGenomesInfo,refDir) // remove definitely, do not include in the workflow
    //chChromSizes = fetch_chrom_sizes(chGenomesInfo,refDir) 
    /*chDACFileRef = downloadDACFile(chGenomesInfo,refDir) 
    
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

    

    chSNPS_ref = downloadSNPRef(chGenomesInfo) 
    //************************************************************************
    //************************************************************************
    

    fastqc(chSampleInfo) 

    //************************************************************************
    //BAM_PROCESSING
    //************************************************************************
    chTrimFiles = trim(chSampleInfo)
    chAlignFiles = align(chTrimFiles,chGenome,chGenomeIndex) 
    chSortedFiles = sort_bam(chAlignFiles)
    lib_complex_preseq(chSortedFiles) //Ver pq está dando erro
    chUniqueFiles = unique_sam(chSortedFiles) 

    chFilteredFiles = quality_filter(chUniqueFiles) 
    chStatsSamtools = createStatsSamtoolsfiltered(chFilteredFiles) 
    chDedupFiles = dedup(chFilteredFiles) 
    //************************************************************************
    //************************************************************************


    chDACFilteredFiles = dac_exclusion(chDedupFiles,chDACFileRef) 
    chIndexFiles = index_sam(chDACFilteredFiles)

    //End Motif and GC content ***********************************************
    chNameSortedFiles = sort_readname_bam(chDACFilteredFiles)
    createMotifGCfile(chNameSortedFiles,chGenome,chGenomeIndex)
    //************************************************************************

    //SNP Fingerprint using SMaSH ************************************************
    chAllIndexFiles = chIndexFiles.collect()
    chAllBAMandBAIIndexFiles = chAllIndexFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }}

    chSMaSHOutout = createSMaSHFingerPrint(chSNPSMaSH,chSNPS_ref,chAllBAMandBAIIndexFiles)
    chSNPSMaSHPlot = createSMaSHFingerPrintPlot(chSMaSHOutout,chSNPSMaSHPyPlot)
    //*****************************************************************************
    
    
    chBedGraphFiles = bam_to_bedgraph(chIndexFiles)
    chBigWig = bedgraph_to_bigwig(chBedGraphFiles,chChromSizes)

    //Pileups ****************************************************************
    chIGVReportsHtml = igv_sample_reports(chBedGraphFiles,chPileUpBED,chGenome,chGenomeIndex).collect()
    chIGVReportMerged = igv_consolidate_report(chIGVReportsHtml,chMultiQCHousekeepingHeader)

    chBigWigAllFiles = chBigWig.collect()
    chBigWigOnlyFiles = chBigWigAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bw') }} // Filter the bw files
    chIGVSession = igv_session(chBigWigOnlyFiles,chIGVFilestoSessions,chGenomesInfo,chPileUpBED)
    //************************************************************************
        
    
    //Fragment Length Distribution *******************************************
    chFragmentsSize = calcFragsLength(chIndexFiles).collect()
    chFragmentsSizeFiles = chFragmentsSize.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.fragment_sizes.txt') }} // Filter the Fragments Size files
    //************************************************************************

    chPeakFiles = peak_bed_graph(chDACFilteredFiles) 
    chBedFiles = bam_to_bed(chDACFilteredFiles) //
    
    chUniqueFrags = unique_frags(chBedFiles).collect() //
    chPeakAllFiles = chPeakFiles.collect() //

    
    //Fragments and peaks Plot *******************************************************
    chNarrowPeakFiles = chPeakAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.narrowPeak') }} // Filter the narrowPeak files
    chFragAndPeaksFilesReport = frags_and_peaks(chNarrowPeakFiles,chUniqueFrags,\
    chMultiQCFragsHeader,chMultiQCPeaksHeader,chReportFragPeaks,chSampleInfo)
    //*********************************************************************************

    //ENRICHMENT *********************************************************************
    chEnrichmentFilesCSV = enrichment(chDACFilteredFiles,chEnrichmentScript).collect()
    chEnrichmentFilesReport = enrichmentReport(chSampleInfo,chEnrichmentFilesCSV,chReportEnrichment).collect()
    chMergedEnrichmentReport = merge_enrichment_reports(chEnrichmentFilesReport,chMultiQCEnrichmentHeader,chMergeReportEnrichment,chSampleInfo).collect()
    

    //Final Report
    chAllPreviousFiles = Channel.fromPath("${workflow.projectDir}/${params.outputFolder}/")
    chFinalReport = multiqc(chIGVReportMerged,chFragmentsSizeFiles,
        chSNPSMaSHPlot,chEnrichmentFilesReport,chFragAndPeaksFilesReport,chMultiQCConfig,chAllPreviousFiles)

    moveSoftFiles(chFinalReport)*/
    
}

