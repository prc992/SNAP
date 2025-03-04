nextflow.enable.dsl=2

include {trim} from '../../modules/trim'
include {align} from './../modules/align'
include {sort_bam} from '../../modules/sort_bam'
include {lib_complex_preseq} from '../../modules/lib_complex_preseq'
include {unique_sam} from '../../modules/unique_sam'
include {quality_filter} from '../../modules/quality_filter'
include {createStatsSamtoolsfiltered} from '../../modules/createStatsSamtoolsfiltered'
include {dedup} from '../../modules/dedup'

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