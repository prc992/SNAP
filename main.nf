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

process multiqc {
    label 'low_cpu_low_mem'
    container = 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'
    publishDir "$path_sample_multiqc", mode : 'copy'
    tag "All Samples" 
    
    input:
    val(_)
    tuple path ("frag_len_hist.txt"),path ("frag_len_mqc.yml")
    path (chFootPrintPDF)
    path (chFragAndPeaks)
    path (chEnrichmentFiles)
    path (configFile)
    tuple val(sampleId), val(enrichment_mark),path(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    file "multiqc_report.html"
    file "multiqc_data/*"

    script:
    """
    multiqc . 
    """
}


process downloadSNPRef {
    label 'low_cpu_low_mem'
    tag "Dowloading - $genome" 
    publishDir "$path_sample_multiqc", mode : 'copy'

    container = "quay.io/biocontainers/wget:1.21.4"
    
    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    input:
    tuple val(genome), val(faGZFile), val(geneAnnotation), val(dacList), val(snp)
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    output:
    file "SNPs.1e5.${genome}.txt"

    script:
    def snpFile = "SNPs.1e5.${genome}.txt"
    
    """
    wget -O ${snpFile} ${snp}
    """
}

process downloadDACFile {

    label 'low_cpu_low_mem'
    tag "Dowloading DAC File - $genome" 
    //publishDir "$refDir", mode : 'copy'

    container = "quay.io/biocontainers/wget:1.21.4"
    
    input:
    tuple val(genome), val(faGZFile), val(geneAnnotation), val(dacList), val(snp)
    path refDir

    output:
    file "${genome}.DAC.bed"

    script:
    def dacFile = "${genome}.DAC.bed"
    def dacFilegz = "${genome}.DAC.bed.gz"
    
    """
    if [ ! -f ${refDir}/${dacFile} ]; then
        wget -O ${refDir}/${dacFilegz} ${dacList}
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
    tuple val(genome), val(faGZFile), val(geneAnnotation), val(dacList), val(snp)
    path refDir

    output:
    file "${genome}.GeneAnotation.gtf"

    script:
    def gtfFile = "${genome}.GeneAnotation.gtf"
    def gtfFilegz = "${genome}.GeneAnotation.gtf.gz"
    

    """
    if [ ! -f ${refDir}/${gtfFile} ]; then
        wget -O ${refDir}/${gtfFilegz} ${geneAnnotation}
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
    tuple val(genome), val(faGZFile), val(geneAnnotation), val(dacList), val(snp)
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
    val(_)

    exec:
    path_sample_align = path_analysis + "/align/" + sampleId

    output:
    path ('*.stats')
    path ('*.idxstats')
    path ('*.flagstat')

    script:
    """
    # Generate stats file
    samtools stats $sampleBam > ${sampleId}.notFiltered.stats

    # Generate idxstats file
    samtools idxstats $sampleBam > ${sampleId}.notFiltered.idxstats

    # Generate flagstat file
    samtools flagstat $sampleBam > ${sampleId}.notFiltered.flagstat
    """
}

process createStatsSamtoolsfiltered {
    label 'low_cpu_low_mem'
    container = 'quay.io/biocontainers/samtools:1.15.1--h1170115_0'
    publishDir "$path_sample_align", mode : 'copy'
    
    tag "Sample - $sampleId" 

    input:
    tuple val(sampleId),val(path_analysis),path(sampleBam)
    val(_)

    exec:
    path_sample_align = path_analysis + "/align/" + sampleId

    output:
    path ('*.stats')
    path ('*.idxstats')
    path ('*.flagstat')
    path ("samtools_stats_mqc_versions.yml")

    script:
    """
    # Generate stats file
    samtools stats $sampleBam > ${sampleId}.AfterFilter.stats

    # Generate idxstats file
    samtools idxstats $sampleBam > ${sampleId}.AfterFilter.idxstats

    # Generate flagstat file
    samtools flagstat $sampleBam > ${sampleId}.AfterFilter.flagstat

    cat <<-END_VERSIONS > samtools_stats_mqc_versions.yml
    "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process quality_filter {
    label 'low_cpu_low_mem'
    container = 'quay.io/biocontainers/samtools:1.15.1--h1170115_0'
    publishDir "$path_sample_align", mode : 'copy'
    
    tag "Sample - $sampleId" 

    input:
    tuple val(sampleId),val(path_analysis),path(sampleBam)
    val(_)

    exec:
    String strBam = sampleId + '.filtered.unique.sorted.bam'
    path_sample_align = path_analysis + "/align/" + sampleId

    output:
    tuple val(sampleId),val(path_analysis),path('*.bam')
    tuple val(sampleId),path ("samtools_QualityFilter_mqc_versions.yml")

    script:
    """
    samtools view -bh -f 3 -F 3844 -q 30 --threads $task.cpus $sampleBam > $strBam

    cat <<-END_VERSIONS > samtools_QualityFilter_mqc_versions.yml
    "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
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
  val(_)

  output:
  path("*.lc_extrap.txt")
  tuple val(sampleId),path ("preseq_mqc_versions.yml")


  exec:
  path_sample_align = path_analysis + "/align/" + sampleId

  script:
  """
  preseq lc_extrap -B $sortedBam > ${sampleId}.lc_extrap.txt

  cat <<-END_VERSIONS > preseq_mqc_versions.yml
  "${task.process}":
    preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
  END_VERSIONS
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
  val(_)

  output:
  path("*fragment_sizes.txt")
  tuple val(sampleId),path ("bamPEFragmentSize_mqc_versions.yml")

  exec:
  path_sample_frags = path_analysis + "/frag/" + sampleId

  script:
  """
  bamPEFragmentSize -b $sortedBam --outRawFragmentLengths ${sampleId}.fragment_sizes.txt

  cat <<-END_VERSIONS > bamPEFragmentSize_mqc_versions.yml
    "${task.process}":
    deeptools: \$(bamPEFragmentSize --version | sed -e "s/bamPEFragmentSize //g")
  END_VERSIONS
  """
}

process fragLenHist {
    container = 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'
    label 'med_cpu_med_mem'
    tag "All Samples"

    publishDir "$path_sample_frags", mode : 'copy'

    input:
    path raw_fragments
    val (_)
    each path (frag_len_header_multiqc)
    each path (chCalcFragHist)
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    exec:
    path_sample_frags = path_analysis + "/frag/"

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
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

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
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

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
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    exec:
    path_sample_multiqc =  path_analysis + "/reports/multiqc/" 

    output:
    path ("*_mqc.csv")

    script:
    """
    python $chMergeReportEnrichment
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

process downloadGenome {

    label 'low_cpu_low_mem'
    tag "Dowloading - $genome" 
    //publishDir "$genomeOut", mode : 'copy'

    container = "quay.io/biocontainers/wget:1.21.4"
    
    exec:
    genomeOut = refDir

    input:
    tuple val(genome), val(faGZFile), val(geneAnnotation), val(dacList), val(snp)
    path refDir

    
    output:
    file "${genome}.fa"

    
    script:
    def genomeFile = "${genome}.fa"
    def genomeFilegz = "${genome}.fa.gz"
    
    """
    if [ ! -f ${refDir}/${genomeFile} ]; then
        wget -O ${refDir}/${genomeFilegz} ${faGZFile}
        gunzip ${refDir}/${genomeFilegz} 
    else
        echo "File ${refDir}/${genomeFile} already exists. Skipping download."
    fi
    ln -s ${refDir}/${genomeFile} ${genomeFile}
    """
}



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
    //chSNPS_ref = Channel.fromPath("$params.snps_ref")
    chMultiQCConfig = Channel.fromPath("$params.multiqc_config")
    chMultiQCFragLenHeader = Channel.fromPath("$params.multiqc_frag_len_header")
    chMultiQCFragPeaksHeader = Channel.fromPath("$params.multiqc_tot_frag_peaks_header")
    chMultiQCEnrichmentHeader = Channel.fromPath("$params.multiqc_enrichment_header")
    
    // Create the genome directory if it doesn't exist
    """
    mkdir -p ${projectDir}/ref_files/genome
    """.execute().waitFor()

    refDir = Channel.fromPath("${projectDir}/ref_files/genome")

    chGenomesSheet = Channel.fromPath(params.genomeInfoPaths)

    chGenomesInfo = chGenomesSheet \
        | splitCsv(header:true) \
        | filter { row -> row.Genome == params.genome } \
        | ifEmpty { error "No matching Genome found in the GenomePaths spreadsheet. Exiting workflow." }
        | map { row-> tuple(row.Genome,row.faGZFile,row.GeneAnotation, row.DACList,row.SNP) }

    // Destructure and store each column into separate variables
    chGenomesInfo
        .map { genome, faGZFile, geneAnnotation, dacList, snp ->
            [genome, faGZFile, geneAnnotation, dacList, snp]
        }


    chGenome = downloadGenome(chGenomesInfo,refDir)
    chGenomeIndex = createGenomeIndex(chGenomesInfo,chGenome,refDir)
    chGeneAnotation = downloadGeneAnotation(chGenomesInfo,refDir)
    chChromSizes = fetch_chrom_sizes(chGenomesInfo,refDir)
    chDACFileRef = downloadDACFile(chGenomesInfo,refDir)
    
    
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


    chSNPS_ref = downloadSNPRef(chGenomesInfo,chSampleInfo)
    fastqc(chSampleInfo) // yaml ready
    chTrimFiles = trim(chSampleInfo) // yaml ready
    chAlignFiles = align(chTrimFiles,chGenome,chGenomeIndex) // yaml ready
    chSortedFiles = sort_bam(chAlignFiles) // yaml ready
    lib_complex(chSortedFiles) // yaml ready
    lib_complex_preseq(chSortedFiles) // yaml ready
    chUniqueFiles = unique_sam(chSortedFiles) // yaml ready
    chStatsSamtools = createStatsSamtools(chUniqueFiles) // yaml ready
    chFilteredFiles = quality_filter(chUniqueFiles) // yaml ready
    chStatsSamtools = createStatsSamtoolsfiltered(chFilteredFiles) // yaml ready
    chDedupFiles = dedup(chFilteredFiles) // yaml ready
    chDACFilteredFiles = dac_exclusion(chDedupFiles,chDACFileRef) // yaml ready

    chIndexFiles = index_sam(chDACFilteredFiles) // yaml ready
    chFragmentsSize = calcFragsLength(chIndexFiles).collect() // yaml ready

    //Verificar se é necessário pois o deepTools já faz isso
    chfragHist = fragLenHist(chFragmentsSize,chMultiQCFragLenHeader,chReportFragHist,chSampleInfo)
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
    chFragAndPeaksFilesReport = frags_and_peaks(chNarrowPeakFiles,chUniqueFrags,chMultiQCFragPeaksHeader,chReportFragPeaks,chSampleInfo)
    //****************************************************************************

    // Processo de SNP Fingerprint
    chSnpFingerprintComplete = snp_fingerprint(chIndexFiles, chSNPS_ref, chGenome).collect()
    chFootPrintPDF = snp_footprint_clustering(chSnpFingerprintComplete,chRSNPFootprint,chSampleInfo)

    
    //ENRICHMENT      ***************************************************
    chEnrichmentFilesCSV = enrichment(chDACFilteredFiles,chEnrichmentScript).collect()
    chEnrichmentFilesReport = enrichmentReport(chSampleInfo,chEnrichmentFilesCSV,chReportEnrichment,chSampleInfo).collect()
    //chEnrichmentFilesReport.subscribe { collectedFiles ->println "Arquivos coletados enrichmentReport: $collectedFiles"}
    chMergedEnrichmentReport = merge_enrichment_reports(chEnrichmentFilesReport,chMultiQCEnrichmentHeader,chMergeReportEnrichment,chSampleInfo).collect()
    //chMergedEnrichmentReport.subscribe { collectedFiles ->println "Arquivos coletados MergedEnrichment: $collectedFiles"}

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

    multiqc(chSnpFingerprintComplete,chfragHist,\
        chFootPrintPDF,chEnrichmentFilesReport,chFragAndPeaksFilesReport,chMultiQCConfig,chSampleInfo)
    
    // COLOCANDO COMO COMENTÁRIO POIS ESTÁ DANDO ERRO POR FALTA DE CONEXÃO
    //pileups_report(chBWFiles,chChromSizes,chPileUpBED,chRPileups)*/


    /*//Collect all files output and the pass to me program that will merge then
    //chAllFiles = chBWFiles.collectFile()
    //pileups_report_comp(chSampleDirPileUps,chChromSizes,chAllFiles,chPileUpBED,chRComparison)*/
}

