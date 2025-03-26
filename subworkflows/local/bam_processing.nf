nextflow.enable.dsl=2

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
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'

process filter_properly_paired {
  label 'low_cpu_high_mem'

  container = params.containers.samtools

  tag "Sample - $sampleId" 
  
  publishDir "${workflow.projectDir}/${params.outputFolder}/align/${sampleId}", mode : 'copy'
  
  input:
  tuple val(sampleId), val(enrichment_mark), val(control), val(read_method), path(sampleBam), val(_)
  
  output:
  tuple val(sampleId), val(enrichment_mark), val(control), val(read_method), path('*.bam'), path("samtools_filter_pp_mqc_versions.yml")

  script:
  String strPPBam = sampleId + '.pp.sorted.bam'

  def filterCommand = ""

  if (read_method == "PE") {
    filterCommand = "samtools view -b -f 2 $sampleBam > ${strPPBam}"
  } else {
    // Just copy the original BAM
    filterCommand = "cp $sampleBam ${strPPBam}"
  }

  """
  echo "Filtering BAM for sample $sampleId (mode: $read_method)"
  $filterCommand

  cat <<-END_VERSIONS > samtools_filter_pp_mqc_versions.yml
  "${task.process}":
      samtools: \$(samtools --version | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}

workflow BAM_PROCESSING {

    take:
    chSampleInfo
    chGenome
    chGenomeIndex
    chChromSizes
    chDACFileRef
    chSNPS_ref
    chAlign
    chSNPSMaSH
    chSNPSMaSHPyPlot
    chMultiQCConfig
    chFilesReportInitialization
    chInitReport
    chFilesReportAlignment
    chAlignmentReport

 
    main:

    if (params.deduped_bam) {
        chDedup = chAlign.map { sampleId,control ,bam, alignYml -> 
            tuple(sampleId, control,bam, null, alignYml)}

        chSortBam = Channel.of("NO_DATA")
        chFilterPP = Channel.of("NO_DATA")
        chLibComplexPreseq = Channel.of("NO_DATA")
        chUniqueSam = Channel.of("NO_DATA")
        chFilteredFiles = Channel.of("NO_DATA")
    }
    else{
        chSortBam = sort_bam(chAlign)
        chFilterPP = filter_properly_paired(chSortBam)
        chLibComplexPreseq = lib_complex_preseq(chFilterPP)
        chUniqueSam = unique_sam(chFilterPP)
        chFilterQuality = quality_filter(chUniqueSam,chSampleInfo)
        chDedup = dedup(chFilterQuality)
    }

    // Filter the DAC files
    if (params.exclude_dac_regions) {
        chDACFilteredFiles = dac_exclusion(chDedup,chDACFileRef)
        //chDACFilteredFiles = chDedup
    } else {
        chDACFilteredFiles = chDedup
    }

    chCreateStatsSamtoolsfiltered = createStatsSamtoolsfiltered(chDACFilteredFiles)
    chIndexFiles = index_sam(chDACFilteredFiles)

    //SNP Fingerprint using SMaSH ************************************************
    chAllIndexFiles = chIndexFiles.collect()
    chAllBAMandBAIIndexFiles = chAllIndexFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }}

    chSMaSHOutout = createSMaSHFingerPrint(chSNPSMaSH,chSNPS_ref,chAllBAMandBAIIndexFiles)
    chSNPSMaSHPlot = createSMaSHFingerPrintPlot(chSMaSHOutout,chSNPSMaSHPyPlot)

    // Collect all the files to generate the MultiQC report
    chSortBamAll = chSortBam.collect()
    chFilterPPAll = chFilterPP.collect()
    chLibComplexPreseqAll = chLibComplexPreseq.collect()
    chUniqueSamAll = chUniqueSam.collect()
    chFilterQualityAll = chFilterQuality.collect()
    chCreateStatsSamtoolsfilteredAll = chCreateStatsSamtoolsfiltered.collect()
    chDedupAll = chDedup.collect()
    chDACFilteredFilesAll = chDACFilteredFiles.collect()
    chSMaSHOutoutAll = chSMaSHOutout.collect()
    chSNPSMaSHPlotAll = chSNPSMaSHPlot.collect()
    chFilesReportInitializationAll = chFilesReportInitialization.collect()
    chFilesReportAlignmentAll = chFilesReportAlignment.collect()

    // Combine all the channels
    chAllChannels = chSortBamAll
        .combine(chSortBamAll)
        .combine(chFilterPPAll)
        .combine(chLibComplexPreseqAll)
        .combine(chUniqueSamAll)
        .combine(chCreateStatsSamtoolsfilteredAll)
        .combine(chFilterQualityAll)
        .combine(chCreateStatsSamtoolsfilteredAll)
        .combine(chDedupAll)
        .combine(chDACFilteredFilesAll)
        .combine(chSMaSHOutoutAll)
        .combine(chSNPSMaSHPlotAll)
        .combine(chFilesReportInitializationAll)
        .combine(chFilesReportAlignmentAll)
    
    // Filter only the files that will be used in the MultiQC report and remove duplicates
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
    chFilesReportBamProcessing = chOnlyFiles.collect()

    // Create the MultiQC report and move the soft files only if this is the last process
    if (params.until == 'BAM_PROCESSING') {
        chBAMProcessReport = multiqc(chInitReport,chFilesReportBamProcessing,chMultiQCConfig)
        moveSoftFiles(chBAMProcessReport)
    } else {
        chBAMProcessReport = Channel.of("NO_DATA")
    }

    emit: bam_processed = chDACFilteredFiles
    emit: bam_processed_index = chIndexFiles
    emit: report_SNP_SMaSH = chSNPSMaSHPlot
    emit: lib_complex = chLibComplexPreseq
    emit: files_report_bam_processing = chFilesReportBamProcessing
    emit: bam_process_report = chBAMProcessReport
}