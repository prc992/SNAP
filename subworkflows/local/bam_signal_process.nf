nextflow.enable.dsl=2

include {call_peaks} from '../../modules/local/call_peaks'
include {peaks_report} from '../../modules/local/peaks_report.nf'
include {bam_to_bedgraph} from '../../modules/local/bam_to_bedgraph'
include {bedgraph_to_bigwig} from '../../modules/local/bedgraph_to_bigwig'
include {igv_reports} from '../../modules/local/igv_reports'
include {igv_sample_reports} from '../../modules/local/igv_reports'
include {igv_consolidate_report} from '../../modules/local/igv_reports'
include {igv_session} from '../../modules/local/igv_reports'
include {enrichmentReport} from '../../modules/local/enrichmentReport'
include {merge_enrichment_reports} from '../../modules/local/merge_enrichment_reports'

workflow BAM_SIGNAL_PROCESSING {

    take:
    chSampleInfo
    chDACFilteredFiles
    chIndexFiles
    chChromSizes
    chPileUpBED
    chGenome
    chGenomeIndex
    chMultiQCHousekeepingHeader
    chIGVFilestoSessions
    chGenomesInfo
    chMultiQCPeaksHeader
    chReportPeaks
    chEnrichmentScript
    chReportEnrichment
    chMergeReportEnrichment
    chMultiQCEnrichmentHeader

    main:
    
    chBedGraphFiles = bam_to_bedgraph(chIndexFiles)
    chBigWig = bedgraph_to_bigwig(chBedGraphFiles,chChromSizes)

    //Pileups ****************************************************************
    chIGVReportsHtml = igv_sample_reports(chBedGraphFiles,chPileUpBED,chGenome,chGenomeIndex).collect()
    chIGVReportMerged = igv_consolidate_report(chIGVReportsHtml,chMultiQCHousekeepingHeader)

    chBigWigAllFiles = chBigWig.collect()
    chBigWigOnlyFiles = chBigWigAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bw') }} // Filter the bw files
    chIGVSession = igv_session(chBigWigOnlyFiles,chIGVFilestoSessions,chGenomesInfo,chPileUpBED)

    chPeakFiles = call_peaks(chDACFilteredFiles) 
    chPeakAllFiles = chPeakFiles.collect()
    chNarrowPeakFiles = chPeakAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.narrowPeak') }} // Filter the narrowPeak files

    chPeaksFilesReport = peaks_report(chNarrowPeakFiles,chMultiQCPeaksHeader,chReportPeaks)

    //ENRICHMENT *********************************************************************
    chEnrichmentFilesCSV = enrichment(chDACFilteredFiles,chEnrichmentScript).collect()
    chEnrichmentFilesReport = enrichmentReport(chSampleInfo,chEnrichmentFilesCSV,chReportEnrichment).collect()
    chMergedEnrichmentReport = merge_enrichment_reports(chEnrichmentFilesReport,chMultiQCEnrichmentHeader,chMergeReportEnrichment,chSampleInfo).collect()

    

}