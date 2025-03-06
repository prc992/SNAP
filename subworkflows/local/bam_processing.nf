nextflow.enable.dsl=2

include {trim} from '../../modules/local/trim'
include {align} from '../../modules/local/align'
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
include {multiqc_bam_processing} from '../../modules/local/multiqc'


workflow BAM_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chChromSizes
    chDACFileRef
    chSNPSMaSH
    chSNPS_ref
    chSNPSMaSHPyPlot
    chMultiQCConfig

    main:

    chTrim = trim(chSampleInfo)
    chAlign = align(chTrim,chGenome,chGenomeIndex)
    chSortBam = sort_bam(chAlign)
    chLibComplexPreseq = lib_complex_preseq(chSortBam)
    chUniqueSam = unique_sam(chSortBam)
    chFilteredFiles = quality_filter(chUniqueSam)
    chCreateStatsSamtoolsfiltered = createStatsSamtoolsfiltered(chFilteredFiles)
    chDedup = dedup(chFilteredFiles)

    // Filter the DAC files
    if (params.exclude_dac_regions) {
        chDACFilteredFiles = dac_exclusion(chDedup,chDACFileRef)
    } else {
        chDACFilteredFiles = chDedup
    }

    chIndexFiles = index_sam(chDACFilteredFiles)

    //SNP Fingerprint using SMaSH ************************************************
    chAllIndexFiles = chIndexFiles.collect()
    chAllBAMandBAIIndexFiles = chAllIndexFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }}

    chSMaSHOutout = createSMaSHFingerPrint(chSNPSMaSH,chSNPS_ref,chAllBAMandBAIIndexFiles)
    chSNPSMaSHPlot = createSMaSHFingerPrintPlot(chSMaSHOutout,chSNPSMaSHPyPlot)


    chTrimAll = chTrim.collect()
    chAlignAll = chAlign.collect()
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

    chAllChannels = chTrimAll
        .combine(chAlignAll)
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
    
    chOnlyFiles = chAllChannels
        .map { values -> values.findAll { it instanceof Path } }
        .flatten() // Garante que os arquivos estejam em um único canal
        .distinct() // Remove arquivos duplicados
    
    multiqc_bam_processing(chOnlyFiles,chMultiQCConfig)

        
    emit: bam_processed = chDACFilteredFiles
    emit: bam_processed_index = chIndexFiles
    emit: report_SNP_SMaSH = chSNPSMaSHPlot
    emit: lib_complex = chLibComplexPreseq
}