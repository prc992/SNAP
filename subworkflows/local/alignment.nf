nextflow.enable.dsl=2

// Import the required processes from the modules
include {trim} from '../../modules/local/trim'
include {align} from '../../modules/local/align'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'

workflow ALIGNMENT {


    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chFilesReportInitialization
    chInitReport
    chMultiQCConfig

    main:

    chTrim = trim(chSampleInfo)
    chAlign = align(chTrim,chGenome,chGenomeIndex)

    // Collect all the files to generate the MultiQC report
    chTrimAll = chTrim.collect()
    chAlignAll = chAlign.collect()

    // Combine all the channels
    chAllChannels = chTrimAll
        .combine(chAlignAll)
        .combine(chFilesReportInitialization)
    
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
    chFilesReportAlignment = chOnlyFiles.collect()

    chAlignmentReport = multiqc(chAlignAll,chFilesReportAlignment,chMultiQCConfig)
    moveSoftFiles(chAlignmentReport)

    emit: align = chAlign
    emit: files_report_alignment = chFilesReportAlignment
    emit: aligment_report = chAlignmentReport

}