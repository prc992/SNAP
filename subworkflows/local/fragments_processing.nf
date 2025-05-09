nextflow.enable.dsl=2

include {sort_readname_bam} from '../../modules/local/sort_bam'
include {bam_to_bed} from '../../modules/local/bam_to_bed'
include {unique_frags} from '../../modules/local/unique_frags'
include {calcFragsLengthDistribuition} from '../../modules/local/calcFragsLength'
include {createMotifGCfile} from '../../modules/local/end_motif_gc'
include {frags_report} from '../../modules/local/frags_report.nf'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'
include {fragle_ct_estimation} from '../../modules/local/fragle_ct_estimation'

process ct_report {
    label 'low_cpu_low_mem'
    container = params.containers.python
    tag "All Samples" 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/multiqc/", mode : 'copy'
    
    input:
    path (chNarrowPeakFiles)
    each path (chMultiQCCTHeader)
    each path (chReportCT)

    output:
    path ("*_mqc.csv")

    script:
    """
    python $chReportCT
    """
}

process filter_bam_fragle {

    label 'low_cpu_low_mem'
    container = params.containers.samtools
    tag "$sampleId"

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/filter_fragle/", mode: 'copy'

    input:
    tuple val(sampleId), val(enrichment_mark), val(control), val(read_method), path(sortedBam), path(sampleBamIndex), val(_)
    each path(chFragleSites)

    output:
    tuple val(sampleId), val(enrichment_mark), val(control), val(read_method), path("*.bam"), path("*.bai"), val(_)

    script:
    """
    SITES_BED="$chFragleSites/$enrichment_mark/sites.bed"

    if [ -f "\$SITES_BED" ]; then
        echo "Filtrando $sortedBam com \$SITES_BED"

        samtools view -b -L "\$SITES_BED" "$sortedBam" > "${sampleId}.filtered.fragle.bam"
        samtools index "${sampleId}.filtered.fragle.bam"
    else
        echo "Arquivo \$SITES_BED não encontrado. Usando arquivos BAM originais."
        cp "$sortedBam" "${sampleId}.filtered.fragle.bam"
        cp "$sampleBamIndex" "${sampleId}.filtered.fragle.bam.bai"
    fi
    """
}

workflow FRAGMENTS_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chBAMProcessedFiles
    chBAMProcessedIndexFiles
    chBAMBAIProcessedFiles
    chMultiQCFragsHeader
    chMultiQCCTHeader
    chReportFrags
    chMultiQCConfig
    chFilesReportInitialization
    chFilesReportBamProcessing
    chBAMProcessReport
    chReportCT

    main:

    
    //End Motif and GC content ***********************************************
    chNameSortedFiles = sort_readname_bam(chBAMProcessedFiles)
    
    chMotifGCfile = createMotifGCfile(chNameSortedFiles,chGenome,chGenomeIndex)
    //************************************************************************
        
    //Fragment Length Distribution *******************************************
    chFragmentsSize = calcFragsLengthDistribuition(chBAMProcessedIndexFiles).collect()
    chFragmentsSizeFiles = chFragmentsSize.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.fragment_sizes.txt') }} // Filter the Fragments Size files
    //************************************************************************

    // Fragle CT estimation **************************************************
    // #######################################################################
    chFragleSites = Channel.fromPath("$params.fragle_sites_ref")

    chAllBAMProcessedIndexFilteredSitesFiles = filter_bam_fragle(chBAMProcessedIndexFiles,chFragleSites).collect()
    chFragleBAMandBAIIndexFiles = chAllBAMProcessedIndexFilteredSitesFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }}
    
    chFragleFiles = fragle_ct_estimation(chFragleBAMandBAIIndexFiles)
    chCTFragleFilesReport = ct_report(chFragleFiles,chMultiQCCTHeader,chReportCT)

    
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
    chCTFragleFilesReportAll = chCTFragleFilesReport.collect()

    // Combine all the channels
    chAllChannels = chNameSortedFilesAll
        .combine(chMotifGCfileAll)
        .combine(chFragmentsSizeFiles)
        .combine(chBedFilesAll)
        .combine(chUniqueFrags)
        .combine(chFragFilesReportAll)
        .combine(chCTFragleFilesReportAll)
        .combine(chFilesReportBamProcessing)
        .combine(chFilesReportInitialization)
        .combine(chFragleFiles)
    
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