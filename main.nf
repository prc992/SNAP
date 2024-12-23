 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

include {lib_complex} from './modules/lib_complex'
include {fastqc} from './modules/fastqc'
include {align} from './modules/align'
include {sort_bam} from './modules/sort_bam'
include {unique_sam} from './modules/unique_sam'
include {enrichment} from './modules/enrichment'
include {index_sam} from './modules/index_sam'
include {dedup} from './modules/dedup'
include {dac_exclusion} from './modules/dac_exclusion'
include {fetch_chrom_sizes} from './modules/fetch_chrom_sizes'
include {peak_bed_graph} from './modules/peak_bed_graph'
include {bam_to_bed} from './modules/bam_to_bed'
include {unique_frags} from './modules/unique_frags'
include {trim} from './modules/trim'
include {snp_fingerprint} from './modules/snp_fingerprint'
include {bedGraphToBigWig} from './modules/bedGraphToBigWig'
include {lenght_fragment_dist_step1} from './modules/lenght_fragment_dist_step'
include {lenght_fragment_dist_step2} from './modules/lenght_fragment_dist_step'
include {pileups_report} from './modules/pileups_report'
include {uropa} from './modules/uropa'
include {snp_footprint_clustering} from './modules/snp_footprint_clustering'

process multiqc_v2 {
    label 'low_cpu_low_mem'
    container = 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'
    publishDir "$path_sample_multiqc", mode : 'copy'
    
    input:
    val(_)
    tuple path ("frag_len_hist.txt"),path ("frag_len_mqc.yml")
    path (chFootPrintPDF)
    path (chFragAndPeaks)
    path (chEnrichmentFiles)
    path (configFile)
    path (chOutputDir)

    exec:
    path_sample_multiqc =  chOutputDir + "/reports/multiqc/" 

    output:
    file "multiqc_report.html"
    file "multiqc_data/*"

    script:
    """
    multiqc . 
    """
}

process downloadGenome {

    label 'low_cpu_low_mem'
    tag "Dowloading - $genome" 
    //publishDir "$genomeOut", mode : 'copy'

    container = "quay.io/biocontainers/wget:1.21.4"
    
    exec:
    genomeOut = refDir

    input:
    val genome
    path refDir

    output:
    file "${genome}.fa"

    
    script:
    def genomeFile = "${genome}.fa"
    def genomeFilegz = "${genome}.fa.gz"
    
    if (genome == 'hg19') {
        url = params.hg19GenomeDownload
    } else if (genome == 'hg38') {
        url = params.hg38GenomeDownload
    } else {
        error "Invalid genome parameter: ${genome}. Allowed values are: ${params.allowedGenomes.join(', ')}"
    }
    """
    if [ ! -f ${refDir}/${genomeFile} ]; then
        wget -O ${refDir}/${genomeFilegz} ${url}
        gunzip ${refDir}/${genomeFilegz} 
    else
        echo "File ${refDir}/${genomeFile} already exists. Skipping download."
    fi
    ln -s ${refDir}/${genomeFile} ${genomeFile}
    """
}

process downloadDACFile {

    label 'low_cpu_low_mem'
    tag "Dowloading DAC File - $genome" 
    //publishDir "$refDir", mode : 'copy'

    container = "quay.io/biocontainers/wget:1.21.4"
    
    input:
    val genome
    path refDir

    output:
    file "${genome}.DAC.bed"

    script:
    def dacFile = "${genome}.DAC.bed"
    def dacFilegz = "${genome}.DAC.bed.gz"
    
    if (genome == 'hg19') {
        url = params.hg19DACListDownload
    } else if (genome == 'hg38') {
        url = params.hg38DACListDownload
    } else {
        error "Invalid genome parameter: ${genome}. Allowed values are: ${params.allowedGenomes.join(', ')}"
    }
    """
    if [ ! -f ${refDir}/${dacFile} ]; then
        wget -O ${refDir}/${dacFilegz} ${url}
        gunzip ${refDir}/${dacFilegz} 
    else
        echo "File ${refDir}/${dacFile} already exists. Skipping download."
    fi
    ln -s ${refDir}/${dacFile} ${dacFile}
    """
}

process downloadGeneAnotation {

    label 'low_cpu_low_mem'
    tag "Dowloading Gene Anotation File - $genome" 
    //publishDir "$refDir", mode : 'copy'

    container = "quay.io/biocontainers/wget:1.21.4"
    
    input:
    val genome
    path refDir

    output:
    file "${genome}.GeneAnotation.gtf"

    script:
    def gtfFile = "${genome}.GeneAnotation.gtf"
    def gtfFilegz = "${genome}.GeneAnotation.gtf.gz"
    
    if (genome == 'hg19') {
        url = params.hg19GeneAnotationDownload
    } else if (genome == 'hg38') {
        url = params.hg38GeneAnotationDownload
    } else {
        error "Invalid genome parameter: ${genome}. Allowed values are: ${params.allowedGenomes.join(', ')}"
    }
    """
    if [ ! -f ${refDir}/${gtfFile} ]; then
        wget -O ${refDir}/${gtfFilegz} ${url}
        gunzip ${refDir}/${gtfFilegz} 
    else
        echo "File ${refDir}/${gtfFile} already exists. Skipping download."
    fi
    ln -s ${refDir}/${gtfFile} ${gtfFile}
    """
}

process createGenomeIndex {
    label 'high_cpu_high_mem'
    container = 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'
    //publishDir "$refDir", mode : 'copy'

    tag "Creating Index - $genome" 

    input:
    val genome
    path genomeFile
    path refDir

    output:
    path "${genome}.fa.*"

    script:
    def genomeFilePac = "${genome}.fa.pac"
    """
    if [ ! -f ${refDir}/${genomeFilePac} ]; then
        echo "${refDir}/${genomeFilePac} file not found. Creating index files."
        bwa index ${genomeFile}
    else
        echo "Index files already exist in ${refDir}. Skipping index creation."
        echo "Creating symlinks to index files."
    fi
    ln -s ${refDir}/${genome}.fa.* .
    """
}

process createSamplesheet {
    label 'low_cpu_low_mem'
    //container = "quay.io/biocontainers/wget:1.21.4"
    tag "Creating Samplesheet" 

    publishDir "$projectDir/$output_dir", mode : 'copy'

    input:
    val sample_dir
    val output_dir
    val enrichment_mark

    output:
    path "snap-samplesheet-*.csv"

    script:
    """
    now=\$(date +'%Y-%m-%d-%H-%M-%S')
    filename="snap-samplesheet-\$now.csv"
    echo "sampleId,enrichment_mark,path,read1,read2" > \$filename

    for subfolder in \$(find ${sample_dir} -mindepth 1 -maxdepth 1 -type d); do
        sampleId=\$(basename \$subfolder)
        # Exclude hidden directories
        if [[ \$sampleId == .* ]]; then
            continue
        fi
        files=(\$(find \$subfolder -type f \\( -name '*.fq.gz' -o -name '*.fq' -o -name '*.fastq.gz' \\) | sort))
        read1=\$(realpath \${files[0]})
        read2=""
        if [ \${#files[@]} -gt 1 ]; then
            read2=\$(realpath \${files[1]})
        fi
        echo "\$sampleId,${enrichment_mark},${output_dir},\$read1,\$read2" >> \$filename
    done
    """
}

process createStatsSamtools {
    label 'low_cpu_low_mem'
    container = 'quay.io/biocontainers/samtools:1.15.1--h1170115_0'
    publishDir "$path_sample_align", mode : 'copy'
    
    tag "Sample - $sampleId" 

    input:
    tuple val(sampleId),val(path_analysis),path(sampleBam)

    exec:
    path_sample_align = path_analysis + "/align/" + sampleId

    output:
    path ('*.stats')
    path ('*.idxstats')
    path ('*.flagstat')

    script:
    """
    # Generate stats file
    samtools stats $sampleBam > ${sampleId}.stats

    # Generate idxstats file
    samtools idxstats $sampleBam > ${sampleId}.idxstats

    # Generate flagstat file
    samtools flagstat $sampleBam > ${sampleId}.flagstat
    """
}

process lib_complex_preseq {
  label 'med_cpu_high_mem'

  //Docker Image
  container = "quay.io/biocontainers/preseq:2.0.2--gsl1.16_0"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_align", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sortedBam)

  output:
  path("*.lc_extrap.txt")

  exec:
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  preseq lc_extrap -B $sortedBam > ${sampleId}.lc_extrap.txt
  """
}

process calcFragsLength {
  label 'med_cpu_med_mem'

  //Docker Image
  //docker pull mgibio/deeptools:3.5.3
  //container = "quay.io/biocontainers/deeptools:2.2.2--py27_0"
  container = "mgibio/deeptools:3.5.3"

  tag "Sample - $sampleId"  
  publishDir "$path_sample_frags", mode : 'copy'

  input:
  tuple val(sampleId),val(path_analysis),path(sortedBam),path (sampleBamIndex)

  output:
  path("*fragment_sizes.txt")

  exec:
  path_sample_frags = path_analysis + "/frag/" + sampleId

  script:
  """
  bamPEFragmentSize -b $sortedBam --outRawFragmentLengths ${sampleId}.fragment_sizes.txt
  """
}

process fragLenHist {
    container = 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'
    label 'med_cpu_med_mem'
    tag "All Samples"

    publishDir "$path_sample_frags", mode : 'copy'

    input:
    path raw_fragments
    each path (frag_len_header_multiqc)
    each path (chCalcFragHist)
    val (chOutputDir)

    exec:
    path_sample_frags = "$chOutputDir" + "/frag/"

    output:
    tuple path ("frag_len_hist.txt"),path ("frag_len_mqc.yml")

    script:
    """
    python $chCalcFragHist \\
        --frag_path "*fragment_sizes.txt" \\
        --output frag_len_hist.txt

    cat $frag_len_header_multiqc frag_len_hist.txt > frag_len_mqc.yml
    """

}

process frags_and_peaks {
    container = 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'
    label 'low_cpu_low_mem'
    tag "All Samples"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    path (narrowPeakFiles)
    path (chPeakAllFiles)
    each path (chMultiQCFragPeaksHeader)
    each path (chReportFragHist)
    val (chOutputDir)

    exec:
    path_sample_multiqc =  chOutputDir + "/reports/multiqc/" 

    output:
    path ("frags_and_peaks_mqc.csv")

    script:
    """
    python $chReportFragHist
    """
}

process enrichmentReport {
    container = 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'
    label 'low_cpu_low_mem'
    tag "All Samples"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    tuple val(sampleId), val(enrichment_mark),val(path),path(read1), val(read2)
    path(csvFiles)
    each path (chReportEnrichment)
    val (chOutputDir)

    exec:
    path_sample_multiqc =  chOutputDir + "/reports/multiqc/" 

    output:
    path ("*_report.csv")

    script:
    """
    python $chReportEnrichment --mark ${enrichment_mark} --samplename ${sampleId}
    """
}

process merge_enrichment_reports {
    container = 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'
    label 'low_cpu_low_mem'
    tag "All Samples"

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    path (chEnrichmentFilesReport)
    each path (chMultiQCEnrichmentHeader)
    each path (chMergeReportEnrichment)
    val (chOutputDir)

    exec:
    path_sample_multiqc =  chOutputDir + "/reports/multiqc/" 

    output:
    path ("*_mqc.csv")

    script:
    """
    python $chMultiQCEnrichmentHeader
    """
}

// RETIRAR ##########################
/*
process createBEDRandomFilesMultiqc{
    container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'
    label 'low_cpu_low_mem'
    tag "$genome"

    publishDir "$path_multiqc", mode : 'copy'

    input:
    val genome
    path chromSizesFile
    val chOutputDir

    output:
    path ("*random.regions.bed")

    exec:
    path_multiqc =  chOutputDir + "/reports/multiqc/" 
    nameFile = genome + ".multiqc.random.regions.bed"

    script:
    """
    bedtools random -g $chromSizesFile -seed 42 > $nameFile
    """
}*/

// RETIRAR ##########################
/*
process deeptoolsComputeMatrix{
    container = 'mgibio/deeptools:3.5.3'
    label 'high_cpu_high_mem'
    tag "All Samples"  

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    path (treat_bw)
    path (bedFile)
    path (chOutputDir)

    exec:
    path_sample_multiqc =  chOutputDir + "/reports/multiqc/" 

    output:
    path ("*computeMatrix.mat.gz")

    script:
    """
         computeMatrix scale-regions \\
        --regionsFileName $bedFile \\
        --scoreFileName $treat_bw \\
        --outFileName AllSamples.computeMatrix.mat.gz \\
        --outFileNameMatrix AllSamples.computeMatrix.vals.mat.tab \\
        --numberOfProcessors $task.cpus
    """
}*/

// RETIRAR ##########################
/*
process deeptoolsPlotCorrelation{
    container = 'mgibio/deeptools:3.5.3'
    label 'med_cpu_med_mem'
    tag "All Samples"  

    publishDir "$path_sample_multiqc", mode : 'copy'

    input:
    path (computeMatrix)
    path (chOutputDir)

    exec:
    path_sample_multiqc =  chOutputDir + "/reports/multiqc/" 

    output:
    path ("*plotCorrelation.pdf")
    path ("*plotCorrelation.mat.tab")

    script:
    """
        plotCorrelation \\
        --corData $computeMatrix \\
        --corMethod spearman --whatToPlot heatmap \\
        --plotFile AllSamples.plotCorrelation.pdf \\
        --outFileCorMatrix AllSamples.plotCorrelation.mat.tab \\
        --skipZeros
    """
}*/



workflow {
    // Static information about the pipeline
    def githubPath = "https://github.com/prc992/SNAP"
    def releaseVersion = "v1.0.15 - LOCAL 60 - 2024-12-10"

    // ASCII art for SNAP
    def asciiArt = """
      ███████╗ ███╗   ██╗ █████╗ ██████╗ 
      ██╔════╝ ████╗  ██║██╔══██╗██╔══██╗
      ███████╗ ██╔██╗ ██║███████║██████╔╝
      ╚════██║ ██║╚██╗██║██╔══██║██╔═══╝ 
      ███████╗ ██║ ╚████║██║  ██║██║     
      ╚══════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     
    """

    // Print the introductory message
    println asciiArt
    println "SNAP pipeline running, created by BacaLab. https://bacalab.dana-farber.org/"
    println "SNAP: Streamlined Nextflow Analysis Pipeline for profiling circulating histone modifications identifies tumor epigenomic signatures in cancer plasma."
    println "GitHub repository: ${githubPath}"
    println "Release version: ${releaseVersion}"


    //Auxiliar code
    chEnrichmentScript= Channel.fromPath("$params.pathEnrichmentScript")
    chRfrag_plotFragDist = Channel.fromPath("$params.pathRfrag_plotFragDist")
    chRComparison = Channel.fromPath("$params.pathRComparison")
    chRPileups= Channel.fromPath("$params.pathRPileups")
    chRSNPFootprint = Channel.fromPath("$params.pathSNPFootprint")
    chReportFragHist = Channel.fromPath("$params.pathReportFragHist")
    chReportFragPeaks = Channel.fromPath("$params.pathReportFragPeaks")
    chReportEnrichment = Channel.fromPath("$params.pathReportEnrichment")
    chMergeReportEnrichment = Channel.fromPath("$params.pathMergeReportEnrichment")

    //Assets
    chPileUpBED = Channel.fromPath("$params.genes_pileup_report")
    chSNPS_ref = Channel.fromPath("$params.snps_ref")
    chMultiQCConfig = Channel.fromPath("$params.multiqc_config")
    chMultiQCFragLenHeader = Channel.fromPath("$params.multiqc_frag_len_header")
    chMultiQCFragPeaksHeader = Channel.fromPath("$params.multiqc_tot_frag_peaks_header")
    chMultiQCEnrichmentHeader = Channel.fromPath("$params.multiqc_enrichment_header")
    
    // Create the genome directory if it doesn't exist
    """
    mkdir -p ${projectDir}/ref_files/genome
    """.execute().waitFor()

    refDir = Channel.fromPath("${projectDir}/ref_files/genome")
    chGenome = downloadGenome(params.genome,refDir)
    chGenomeIndex = createGenomeIndex(params.genome,chGenome,refDir)
    chGeneAnotation = downloadGeneAnotation(params.genome,refDir)
    chChromSizes = fetch_chrom_sizes(params.genome,refDir)
    chDACFileRef = downloadDACFile(params.genome,refDir)
    
    


    // If the 'samplesheet' parameter is provided, use it directly; otherwise, create a new samplesheet
    if (params.samplesheet) {
        //println "Using provided samplesheet: ${params.samplesheet}"
        chSampleSheet = Channel.fromPath(params.samplesheet)
    } else {
        //println "Creating samplesheet because none was provided."
        chSampleSheet = createSamplesheet(
            params.sample_dir, 
            params.output_dir, 
            params.enrichment_mark ?: 'no_enrichment_mark'
        )
    }
    
    chSampleInfo = chSampleSheet \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.enrichment_mark,"${projectDir}/${row.path}", row.read1, row.read2) }


    //Extract outputdir from the first row
    chOutputDir = chSampleInfo.first().map { firstItem -> firstItem[2] }

    //RETIRAR
    //chBEDRandomFilesMultiqc = createBEDRandomFilesMultiqc(params.genome,chChromSizes,chOutputDir)
    
    fastqc(chSampleInfo)
    chTrimFiles = trim(chSampleInfo)
    chAlignFiles = align(chTrimFiles,chGenome,chGenomeIndex)    
    chSortedFiles = sort_bam(chAlignFiles)
    lib_complex(chSortedFiles)
    lib_complex_preseq(chSortedFiles)
    chUniqueFiles = unique_sam(chSortedFiles)
    chStatsSamtools = createStatsSamtools(chUniqueFiles)
    chDedupFiles = dedup(chUniqueFiles)
    chDACFilteredFiles = dac_exclusion(chDedupFiles,chDACFileRef)

    chIndexFiles = index_sam(chDACFilteredFiles)
    chFragmentsSize = calcFragsLength(chIndexFiles).collect()

    //Verificar se é necessário pois o deepTools já faz isso
    chfragHist = fragLenHist(chFragmentsSize,chMultiQCFragLenHeader,chReportFragHist,chOutputDir)
    //************************************************************************

    chPeakFiles = peak_bed_graph(chDACFilteredFiles)

    
    uropa(chPeakFiles,chGeneAnotation)
    chBedFiles = bam_to_bed(chDACFilteredFiles)

    chUniqueFrags = unique_frags(chBedFiles).collect()
    chPeakAllFiles = chPeakFiles.collect()

    // Filter the narrowPeak files
    chNarrowPeakFiles = chPeakAllFiles.map { collectedFiles ->
    collectedFiles.findAll { it.toString().endsWith('.narrowPeak') }}
    //In case you want to print the files
    //chNarrowPeakFiles.subscribe { collectedFiles ->println "Arquivos coletados: $collectedFiles"}


    //FRAGMENTS AND PEAKS      ***************************************************
    chFragAndPeaksFilesReport = frags_and_peaks(chNarrowPeakFiles,chUniqueFrags,chMultiQCFragPeaksHeader,chReportFragPeaks,chOutputDir)
    //****************************************************************************

    // Processo de SNP Fingerprint
    chSnpFingerprintComplete = snp_fingerprint(chIndexFiles, chSNPS_ref, chGenome).collect()
    chFootPrintPDF = snp_footprint_clustering(chSnpFingerprintComplete,chRSNPFootprint,chOutputDir)

    
    //ENRICHMENT      ***************************************************
    chEnrichmentFilesCSV = enrichment(chDACFilteredFiles,chEnrichmentScript).collect()
    chEnrichmentFilesReport = enrichmentReport(chSampleInfo,chEnrichmentFilesCSV,chReportEnrichment,chOutputDir).collect()
    chEnrichmentFilesReport.subscribe { collectedFiles ->println "Arquivos coletados enrichmentReport: $collectedFiles"}
    chMergedEnrichmentReport = merge_enrichment_reports(chEnrichmentFilesReport,chMultiQCEnrichmentHeader,chMergeReportEnrichment,chOutputDir).collect()
    chMergedEnrichmentReport.subscribe { collectedFiles ->println "Arquivos coletados MergedEnrichment: $collectedFiles"}
    /*
    //Verificar se é necessário pois o deepTools já faz isso
    chFragDis = lenght_fragment_dist_step1(chDACFilteredFiles)
    lenght_fragment_dist_step2(chFragDis,chRfrag_plotFragDist)
    //************************************************************************

    chBWFiles = bedGraphToBigWig(chPeakFiles,chChromSizes)

    // RETIRAR ##########################
    // DEEPTOOLS_COMPUTEMATRIX
    //chBWAllFiles = chBWFiles.collect()
    //chBWTreatFiles = chBWAllFiles.map { collectedFiles ->
    //collectedFiles.findAll { it.toString().endsWith('treat_pileup.bdg.bw') }}
    //chDeepToolsMatrix = deeptoolsComputeMatrix(chBWTreatFiles,chBEDRandomFilesMultiqc,chOutputDir)

    // DEEPTOOLS_PLOTCORRELATION
    //chPlotCorrelation = deeptoolsPlotCorrelation(chDeepToolsMatrix,chOutputDir)
    // RETIRAR ##########################

    multiqc_v2(chSnpFingerprintComplete,chfragHist,\
        chFootPrintPDF,chEnrichmentFilesReport,chFragAndPeaksFilesReport,chMultiQCConfig,chOutputDir)
    
    // COLOCANDO COMO COMENTÁRIO POIS ESTÁ DANDO ERRO POR FALTA DE CONEXÃO
    pileups_report(chBWFiles,chChromSizes,chPileUpBED,chRPileups)*/


    /*//Collect all files output and the pass to me program that will merge then
    //chAllFiles = chBWFiles.collectFile()
    //pileups_report_comp(chSampleDirPileUps,chChromSizes,chAllFiles,chPileUpBED,chRComparison)*/
}

