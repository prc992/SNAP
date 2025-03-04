nextflow.enable.dsl=2

include {bam_to_bedgraph} from '../../modules/local/bam_to_bedgraph'
include {bedgraph_to_bigwig} from '../../modules/local/bedgraph_to_bigwig'
include {igv_reports} from '../../modules/local/igv_reports'
include {igv_sample_reports} from '../../modules/local/igv_reports'
include {igv_consolidate_report} from '../../modules/local/igv_reports'
include {igv_session} from '../../modules/local/igv_reports'

workflow BAM_SIGNAL_PROCESSING {

    take:
    chDACFilteredFiles
    chIndexFiles
    chChromSizes
    chPileUpBED
    chGenome
    chGenomeIndex
    chMultiQCHousekeepingHeader
    chIGVFilestoSessions
    chGenomesInfo

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

}