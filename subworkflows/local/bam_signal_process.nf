nextflow.enable.dsl=2

include {call_peaks} from '../../modules/local/call_peaks'
include {peaks_report} from '../../modules/local/peaks_report.nf'
include {bam_to_bedgraph} from '../../modules/local/bam_to_bedgraph'
include {bedgraph_to_bigwig} from '../../modules/local/bedgraph_to_bigwig'
include {igv_reports} from '../../modules/local/igv_reports'
include {igv_sample_reports} from '../../modules/local/igv_reports'
include {igv_consolidate_report} from '../../modules/local/igv_reports'
include {igv_session} from '../../modules/local/igv_reports'
include {enrichment} from '../../modules/local/enrichment'
include {enrichmentReport} from '../../modules/local/enrichmentReport'
include {merge_enrichment_reports} from '../../modules/local/merge_enrichment_reports'
include {multiqc_bam_signal_processing} from '../../modules/local/multiqc'

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
    chFilesReportBamProcessing
    chFilesReportInitialization
    chMultiQCConfig

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

    // Collect all the files to generate the MultiQC report
    chBedGraphFilesAll = chBedGraphFiles.collect()
    chBigWigAllFiles = chBigWig.collect()
    chIGVReportsHtmlAll = chIGVReportsHtml.collect()
    chIGVReportMergedAll = chIGVReportMerged.collect()
    chBigWigAllFilesAll = chBigWigAllFiles.collect()
    chIGVSessionAll = chIGVSession.collect()
    chPeakFilesAll = chPeakFiles.collect()
    chPeaksFilesReportAll = chPeaksFilesReport.collect()
    chEnrichmentFilesCSVAll = chEnrichmentFilesCSV.collect()
    chEnrichmentFilesReportAll = chEnrichmentFilesReport.collect()
    chMergedEnrichmentReportAll = chMergedEnrichmentReport.collect()


    // Combine all the channels
    chAllChannelsProcessing = chBedGraphFilesAll
        .combine(chBigWigAllFiles)
        .combine(chIGVReportsHtmlAll)
        .combine(chIGVReportMergedAll)
        .combine(chBigWigAllFilesAll)
        .combine(chIGVSessionAll)
        .combine(chPeakFilesAll)
        .combine(chPeaksFilesReportAll)
        .combine(chEnrichmentFilesCSVAll)
        .combine(chEnrichmentFilesReportAll)
        .combine(chMergedEnrichmentReportAll)
        .combine(chFilesReportBamProcessing)
        .combine(chFilesReportInitialization)
    
    chOnlyFilesProcessing = chAllChannelsProcessing
    .flatten() // Garante que os arquivos estejam em um único fluxo
    .collect() // Junta todos os arquivos antes de processá-los
    .map { files -> 
        def uniqueFiles = [:] as LinkedHashMap
        files.findAll { it instanceof Path } // 🔹 Mantém apenas arquivos (Path)
             .each { file -> uniqueFiles.putIfAbsent(file.getName(), file) } // Mantém apenas a primeira ocorrência do nome
        return uniqueFiles.values() // Retorna apenas os arquivos únicos
    }
    .flatten() // Garante que cada arquivo seja emitido separadamente no canal

    chFilesReportSignalProcess = chOnlyFilesProcessing.collect()
    multiqc_bam_signal_processing(chFilesReportSignalProcess,chMultiQCConfig)



    emit: igv_report_merged = chIGVReportMerged
    emit: merge_enrichment_reports = chMergedEnrichmentReport
    emit: peaks_report = chPeaksFilesReport
    emit: files_report_bam_signal_processing = chFilesReportSignalProcess
}