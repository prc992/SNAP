nextflow.enable.dsl=2

// Import the required processes from the modules
include {trim} from '../../modules/local/trim'
include {trim_fastp} from '../../modules/local/trim'
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
    if (params.trim_method == 'FASTP') {
        println "Using FASTP for trimming"
        chTrim = trim_fastp(chSampleInfo)
    } else {
        println "Using TRIM_GALORE for trimming"
        chTrim = trim(chSampleInfo)
    }

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

    // Create the MultiQC report and move the soft files only if this is the last process
    if (params.until == 'ALIGNMENT') {
        chAlignmentReport = multiqc(chAlignAll,chFilesReportAlignment,chMultiQCConfig)
        moveSoftFiles(chAlignmentReport)
    } else {
        chAlignmentReport = Channel.of("NO_DATA")
    }

    emit: align = chAlign
    emit: files_report_alignment = chFilesReportAlignment
    emit: aligment_report = chAlignmentReport

}