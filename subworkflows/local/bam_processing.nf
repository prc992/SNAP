nextflow.enable.dsl=2

include {trim} from '../../modules/local/trim'
include {align} from '../../modules/local/align'
include {sort_bam} from '../../modules/local/sort_bam'
include {lib_complex_preseq} from '../../modules/local/lib_complex_preseq'
include {unique_sam} from '../../modules/local/unique_sam'
include {quality_filter} from '../../modules/local/quality_filter'
include {createStatsSamtoolsfiltered} from '../../modules/local/createStatsSamtoolsfiltered'
include {dedup} from '../../modules/local/dedup'
include {dac_exclusion} from '../../modules/dac_exclusion'

workflow BAM_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chChromSizes
    chDACFileRef

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

    emit: bam_processed = chDACFilteredFiles
}