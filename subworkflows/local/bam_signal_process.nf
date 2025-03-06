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
    chIGVReportMergedAll = chIGVReportMerged.collect()
    chPeakFilesAll = chPeakFiles.collect()
    chFilesReportBamProcessingAll = chFilesReportBamProcessing.collect()
    chFilesReportInitializationAll = chFilesReportInitialization.collect()

    // Combine all the channels
    chAllChannels = chBedGraphFilesAll
        .combine(chIGVReportMergedAll)
        .combine(chPeakFilesAll)
        .combine(chPeaksFilesReportAll)
        .combine(chIGVReportsHtml)
        .combine(chBigWigAllFiles)
        .combine(chMergedEnrichmentReport)
        .combine(chFilesReportBamProcessingAll)
        .combine(chFilesReportInitializationAll)
    
    // Filter only the files that will be used in the MultiQC report and remove duplicates
    chOnlyFiles = chAllChannels
        .map { values -> 
            values.findAll { 
                it instanceof Path && ( 
                    it.toString().endsWith(".yml") || 
                    it.toString().endsWith(".csv") || 
                    it.toString().endsWith(".zip") || 
                    it.toString().endsWith(".txt") || 
                    it.toString().endsWith(".stats") || 
                    it.toString().endsWith(".txt") || 
                    it.toString().endsWith(".idxstats") ||
                    it.toString().endsWith(".flagstat") ||  
                    it.toString().contains("Dendrogram_of_Samples")
                )
            }
        }
        .flatten() // Garante que os arquivos estejam em um único fluxo
        .reduce( [:] as LinkedHashMap ) { acc, file -> 
            acc.putIfAbsent(file.getName(), file) // Mantém apenas a primeira ocorrência do nome do arquivo
            acc
        }
        .map { it.values().toList() } // 🔹 Converte para uma lista
        chFilesReportSignalProcess = chOnlyFiles.collect()
    
    multiqc_bam_signal_processing(chFilesReportSignalProcess,chMultiQCConfig)



    emit: igv_report_merged = chIGVReportMerged
    emit: merge_enrichment_reports = chMergedEnrichmentReport
    emit: peaks_report = chPeaksFilesReport

    

}