nextflow.enable.dsl=2

include {sort_readname_bam} from '../../modules/local/sort_bam'
include {bam_to_bed} from '../../modules/local/bam_to_bed'
include {unique_frags} from '../../modules/local/unique_frags'
include {calcFragsLengthDistribuition} from '../../modules/local/calcFragsLength'
include {createMotifGCfile} from '../../modules/local/end_motif_gc'
include {frags_report} from '../../modules/local/frags_report.nf'


workflow FRAGMENTS_PROCESSING {

    take:
    chBAMProcessedFiles
    chIndexFiles
    chGenome
    chGenomeIndex
    chMultiQCFragsHeader
    chReportFragPeaks

    main:
    
    //End Motif and GC content ***********************************************
    chNameSortedFiles = sort_readname_bam(chBAMProcessedFiles)
    createMotifGCfile(chNameSortedFiles,chGenome,chGenomeIndex)
    //************************************************************************
        
    //Fragment Length Distribution *******************************************
    chFragmentsSize = calcFragsLengthDistribuition(chIndexFiles).collect()
    chFragmentsSizeFiles = chFragmentsSize.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.fragment_sizes.txt') }} // Filter the Fragments Size files
    //************************************************************************
    
    chBedFiles = bam_to_bed(chBAMProcessedFiles) 

    chUniqueFrags = unique_frags(chBedFiles).collect()
    chFragAndPeaksFilesReport = frags_report(chUniqueFrags,chMultiQCFragsHeader,chReportFragPeaks)
}