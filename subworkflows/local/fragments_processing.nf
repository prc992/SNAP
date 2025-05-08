nextflow.enable.dsl=2

include {sort_readname_bam} from '../../modules/local/sort_bam'
include {bam_to_bed} from '../../modules/local/bam_to_bed'
include {unique_frags} from '../../modules/local/unique_frags'
include {calcFragsLengthDistribuition} from '../../modules/local/calcFragsLength'
include {createMotifGCfile} from '../../modules/local/end_motif_gc'
include {frags_report} from '../../modules/local/frags_report.nf'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'

process fragle_ct_estimation {

    label 'high_cpu_high_plus_mem'
    container = params.containers.fragle

    tag "All Samples"
    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/fragle/", mode : 'copy'
    
    input:
    tuple val(sampleId),val(enrichment_mark),val(control),val(read_method),path(sampleBam),val(_)
    
    output:
    path ("*.csv")
    
    script:
    """
    cd /usr/src/app
    python /usr/src/app/main.py --input \$NXF_WORK --output \$NXF_WORK --mode R --cpu ${task.cpus} --threads ${task.cpus}
    """
}


workflow FRAGMENTS_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chBAMProcessedFiles
    chBAMProcessedIndexFiles
    chMultiQCFragsHeader
    chReportFrags
    chMultiQCConfig
    chFilesReportInitialization
    chFilesReportBamProcessing
    chBAMProcessReport

    main:

    chFragleFiles = fragle_ct_estimation(chBAMProcessedFiles)
    
    //End Motif and GC content ***********************************************
    chNameSortedFiles = sort_readname_bam(chBAMProcessedFiles)
    
    chMotifGCfile = createMotifGCfile(chNameSortedFiles,chGenome,chGenomeIndex)
    //************************************************************************
        
    //Fragment Length Distribution *******************************************
    chFragmentsSize = calcFragsLengthDistribuition(chBAMProcessedIndexFiles).collect()
    chFragmentsSizeFiles = chFragmentsSize.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.fragment_sizes.txt') }} // Filter the Fragments Size files
    //************************************************************************
    
    chBedFiles = bam_to_bed(chNameSortedFiles) 

    chUniqueFrags = unique_frags(chBedFiles).collect()
    chFragFilesReport = frags_report(chUniqueFrags,chMultiQCFragsHeader,chReportFrags)
    
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
        .combine(chFragmentsSizeFiles)
        .combine(chBedFilesAll)
        .combine(chUniqueFrags)
        .combine(chFragFilesReportAll)
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

    chFilesReportFragmentsProcess = chOnlyFilesProcessing.collect()

    // Create the MultiQC report and move the soft files only if this is the last process
    if (params.until == 'FRAGMENTS_PROCESSING') {
        chFragsProcessReport = multiqc(chBAMProcessReport,chFilesReportFragmentsProcess,chMultiQCConfig)
        moveSoftFiles(chFragsProcessReport)
    } else {
        chFragsProcessReport = Channel.of("NO_DATA")
    }

    emit: frag_size_files = chFragmentsSizeFiles
    emit: frag_report = chFragFilesReport
    emit: frag_process_report = chFragsProcessReport
    emit: files_report_fragments_processing = chFilesReportFragmentsProcess
}