nextflow.enable.dsl=2

include {sort_readname_bam} from '../../modules/local/sort_bam'
include {bam_to_bed} from '../../modules/local/bam_to_bed'
include {unique_frags} from '../../modules/local/unique_frags'
include {calcFragsLengthDistribuition} from '../../modules/local/calcFragsLength'
include {createMotifGCfile} from '../../modules/local/end_motif_gc'
include {frags_report} from '../../modules/local/frags_report.nf'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'


workflow FRAGMENTS_PROCESSING {

    take:
    chGenome
    chGenomeIndex
    chBAMProcessedFiles
    chBAMProcessedIndexFiles
    chMultiQCFragsHeader
    chReportFrags
    chMultiQCConfig
    chFilesReportInitialization
    chFilesReportBamProcessing
    chFilesReportSignalProcess
    chBAMSignalReport

    main:
    
    //End Motif and GC content ***********************************************
    chNameSortedFiles = sort_readname_bam(chBAMProcessedFiles)
    chMotifGCfile = createMotifGCfile(chNameSortedFiles,chGenome,chGenomeIndex)
    //************************************************************************
        
    //Fragment Length Distribution *******************************************
    chFragmentsSize = calcFragsLengthDistribuition(chIndexFiles).collect()
    chFragmentsSizeFiles = chFragmentsSize.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.fragment_sizes.txt') }} // Filter the Fragments Size files
    //************************************************************************
    
    chBedFiles = bam_to_bed(chBAMProcessedFiles) 

    chUniqueFrags = unique_frags(chBedFiles).collect()
    chFragFilesReport = frags_report(chUniqueFrags,chMultiQCFragsHeader,chReportFragPeaks)

    // Collect all the files to generate the MultiQC report
    chNameSortedFilesAll = chNameSortedFiles.collect()
    chMotifGCfileAll = chMotifGCfile.collect()
    //chFragmentsSize
    chBedFilesAll = chBedFiles.collect()
    //chUniqueFrags
    chFragFilesReportAll = chFragFilesReport.collect()

    // Combine all the channels
    chAllChannels = chNameSortedFilesAll
        .combine(chMotifGCfileAll)
        .combine(chFragmentsSize)
        .combine(chBedFilesAll)
        .combine(chUniqueFrags)
        .combine(chFragFilesReportAll)
        .combine(chFilesReportSignalProcess)
        .combine(chFilesReportBamProcessing)
        .combine(chFilesReportInitialization)
    
    chOnlyFilesProcessing = chAllChannels
    .flatten() // Garante que os arquivos estejam em um único fluxo
    .collect() // Junta todos os arquivos antes de processá-los
    .map { files -> 
        def uniqueFiles = [:] as LinkedHashMap
        files.findAll { it instanceof Path } // 🔹 Mantém apenas arquivos (Path)
             .each { file -> uniqueFiles.putIfAbsent(file.getName(), file) } // Mantém apenas a primeira ocorrência do nome
        return uniqueFiles.values() // Retorna apenas os arquivos únicos
    }
    .flatten() // Garante que cada arquivo seja emitido separadamente no canal

    chFilesReportFragmentslProcess = chOnlyFilesProcessing.collect()
    chFragsProcessReport = multiqc(chBAMSignalReport,chFilesReportFragmentslProcess,chMultiQCConfig)
    moveSoftFiles(chFragsProcessReport)

    emit: frag_size_files = chFragmentsSizeFiles
    emit: frag_report = chFragFilesReport
    emit: frag_process_report = chFragsProcessReport
}