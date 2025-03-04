nextflow.enable.dsl=2

include {trim} from '../../modules/local/trim'
include {align} from '../../modules/local/align'
include {sort_bam} from '../../modules/local/sort_bam'
include {lib_complex_preseq} from '../../modules/local/lib_complex_preseq'
include {unique_sam} from '../../modules/local/unique_sam'
include {quality_filter} from '../../modules/local/quality_filter'
include {createStatsSamtoolsfiltered} from '../../modules/local/createStatsSamtoolsfiltered'
include {dedup} from '../../modules/local/dedup'

workflow BAM_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex

    main:

    chTrim = trim(chSampleInfo)
    chAlign = align(chTrim, chGenomeIndex)
    chSortBam = sort_bam(chAlign)
    chLibComplexPreseq = lib_complex_preseq(chSortBam)
    chUniqueSam = unique_sam(chLibComplexPreseq)
    chQualityFilter = quality_filter(chUniqueSam)
    chCreateStatsSamtoolsfiltered = createStatsSamtoolsfiltered(chQualityFilter, chGenomeIndex, chChromSizes)
    chDedup = dedup(chCreateStatsSamtoolsfiltered)

    emit: bam_deduped = chDedup
}