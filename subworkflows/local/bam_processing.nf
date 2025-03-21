nextflow.enable.dsl=2

include {sort_bam} from '../../modules/local/sort_bam'
include {lib_complex_preseq} from '../../modules/local/lib_complex_preseq'
include {unique_sam} from '../../modules/local/unique_sam'
include {quality_filter} from '../../modules/local/quality_filter'
include {createStatsSamtoolsfiltered} from '../../modules/local/createStatsSamtoolsfiltered'
include {dedup} from '../../modules/local/dedup'
include {dac_exclusion} from '../../modules/local/dac_exclusion'
include {index_sam} from '../../modules/local/index_sam'
include {createSMaSHFingerPrint} from '../../modules/local/snp_smash_fingerprint'
include {createSMaSHFingerPrintPlot} from '../../modules/local/snp_smash_fingerprint'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'


workflow BAM_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chChromSizes
    chDACFileRef
    chSNPS_ref
    chAlign
    chSNPSMaSH
    chSNPSMaSHPyPlot
    chMultiQCConfig
    chFilesReportInitialization
    chInitReport
    chFilesReportAlignment
    chAlignmentReport

 
    main:

    if (params.deduped_bam) {
        chDedup = chAlign.map { sampleId, bam, alignYml -> 
            tuple(sampleId, bam, null, alignYml)}

        chSortBam = Channel.of("NO_DATA")
        chLibComplexPreseq = Channel.of("NO_DATA")
        chUniqueSam = Channel.of("NO_DATA")
        chFilteredFiles = Channel.of("NO_DATA")
    }
    else{
        println "Arquivos não deduplicados"
        chAlign.view()
        chSortBam = sort_bam(chAlign)
        chLibComplexPreseq = lib_complex_preseq(chSortBam)
        chUniqueSam = unique_sam(chSortBam)
        chFilteredFiles = quality_filter(chUniqueSam)
        chDedup = dedup(chFilteredFiles)
    }

    // Filter the DAC files
    if (params.exclude_dac_regions) {
        chDACFilteredFiles = dac_exclusion(chDedup,chDACFileRef)
        //chDACFilteredFiles = chDedup
    } else {
        chDACFilteredFiles = chDedup
    }

    chCreateStatsSamtoolsfiltered = createStatsSamtoolsfiltered(chDACFilteredFiles)
    chIndexFiles = index_sam(chDACFilteredFiles)

    //SNP Fingerprint using SMaSH ************************************************
    chAllIndexFiles = chIndexFiles.collect()
    chAllBAMandBAIIndexFiles = chAllIndexFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }}

    chSMaSHOutout = createSMaSHFingerPrint(chSNPSMaSH,chSNPS_ref,chAllBAMandBAIIndexFiles)
    chSNPSMaSHPlot = createSMaSHFingerPrintPlot(chSMaSHOutout,chSNPSMaSHPyPlot)

    // Collect all the files to generate the MultiQC report
    chSortBamAll = chSortBam.collect()
    chLibComplexPreseqAll = chLibComplexPreseq.collect()
    chUniqueSamAll = chUniqueSam.collect()
    chFilteredFilesAll = chFilteredFiles.collect()
    chCreateStatsSamtoolsfilteredAll = chCreateStatsSamtoolsfiltered.collect()
    chFilteredFilesAll = chFilteredFiles.collect()
    chCreateStatsSamtoolsfilteredAll = chCreateStatsSamtoolsfiltered.collect()
    chDedupAll = chDedup.collect()
    chDACFilteredFilesAll = chDACFilteredFiles.collect()
    chSMaSHOutoutAll = chSMaSHOutout.collect()
    chSNPSMaSHPlotAll = chSNPSMaSHPlot.collect()
    chFilesReportInitializationAll = chFilesReportInitialization.collect()
    chFilesReportAlignmentAll = chFilesReportAlignment.collect()

    // Combine all the channels
    chAllChannels = chSortBamAll
        .combine(chSortBamAll)
        .combine(chLibComplexPreseqAll)
        .combine(chUniqueSamAll)
        .combine(chFilteredFilesAll)
        .combine(chCreateStatsSamtoolsfilteredAll)
        .combine(chFilteredFilesAll)
        .combine(chCreateStatsSamtoolsfilteredAll)
        .combine(chDedupAll)
        .combine(chDACFilteredFilesAll)
        .combine(chSMaSHOutoutAll)
        .combine(chSNPSMaSHPlotAll)
        .combine(chFilesReportInitializationAll)
        .combine(chFilesReportAlignmentAll)
    
    // Filter only the files that will be used in the MultiQC report and remove duplicates
    chOnlyFiles = chAllChannels
    .flatten() // Make sure the files are in a single flow
    .collect() // Joins all files before processing them
    .map { files -> 
        def uniqueFiles = [:] as LinkedHashMap
        files.findAll { it instanceof Path } // Keeps only files (Path)
             .each { file -> uniqueFiles.putIfAbsent(file.getName(), file) } // Keeps only the first occurrence of the name
        return uniqueFiles.values()  // Returns only unique files
    } 
    .flatten()
    chFilesReportBamProcessing = chOnlyFiles.collect()
    
    chBAMProcessReport = multiqc(chInitReport,chFilesReportBamProcessing,chMultiQCConfig)
    moveSoftFiles(chBAMProcessReport)

    emit: bam_processed = chDACFilteredFiles
    emit: bam_processed_index = chIndexFiles
    emit: report_SNP_SMaSH = chSNPSMaSHPlot
    emit: lib_complex = chLibComplexPreseq
    emit: files_report_bam_processing = chFilesReportBamProcessing
    emit: bam_process_report = chBAMProcessReport
}