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
include {json_uropa} from './modules/uropa'
include {snp_footprint_clustering} from './modules/snp_footprint_clustering'

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

workflow {
    // Static information about the pipeline
    def githubPath = "https://github.com/prc992/SNAP"
    def releaseVersion = "v1.0.3 - LOCAL - 2021-09-30"

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
    println "SNAP: Streamlined Nextflow Analysis Pipeline for profiling circulating histone modifications identifies tumor epigenomic signatures in cancer plasma"
    println "GitHub repository: ${githubPath}"
    println "Release version: ${releaseVersion}"

    chSampleInfo = Channel.fromPath(params.samples) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.path, row.read1, row.read2) }

    //Auxiliar code
    chEnrichmentScript= Channel.fromPath("$params.pathEnrichmentScript")
    chRfrag_plotFragDist = Channel.fromPath("$params.pathRfrag_plotFragDist")
    chRComparison = Channel.fromPath("$params.pathRComparison")
    chRPileups= Channel.fromPath("$params.pathRPileups")
    chRSNPFootprint = Channel.fromPath("$params.pathSNPFootprint")
    //chJson_file = Channel.fromPath("$params.pathJson_file")


    //Assets
    chPileUpBED = Channel.fromPath("$params.genes_pileup_report")
    chSNPS_ref = Channel.fromPath("$params.snps_ref")
    
    // Create the genome directory if it doesn't exist
    """
    mkdir -p ${projectDir}/ref_files/genome
    """.execute().waitFor()

    refDir = Channel.fromPath("${projectDir}/ref_files/genome")
    chGenome = downloadGenome(params.genome,refDir)
    chGenomeIndex = createGenomeIndex(params.genome,chGenome,refDir)

    fastqc(chSampleInfo)
    chTrimFiles = trim(chSampleInfo)
    chAlignFiles = align(chTrimFiles,chSampleInfo,chGenome,chGenomeIndex)
    
    chSortedFiles = sort_bam(chAlignFiles,chSampleInfo)
    lib_complex(chSortedFiles,chSampleInfo)
    chUniqueFiles = unique_sam(chSortedFiles,chSampleInfo)

    chDACFileRef = downloadDACFile(params.genome,refDir)

    
    chDedupFiles = dedup(chUniqueFiles,chSampleInfo)
    chDACFilteredFiles = dac_exclusion(chDedupFiles,chSampleInfo,chDACFileRef)

    chIndexFiles = index_sam(chDACFilteredFiles,chSampleInfo)
    chPeakFiles = peak_bed_graph(chDACFilteredFiles,chSampleInfo)
    /*
    //corrigir depois
    //chJson_file = json_uropa(chSampleInfo)
    //uropa(chPeakFiles,chJson_file,chGTF_ref,chSampleInfo)

    chBedFiles = bam_to_bed(chDedupFiles,chSampleInfo)
    unique_frags(chBedFiles,chSampleInfo)
    chChromSizes = fetch_chrom_sizes(chSampleInfo)
    //snp_fingerprint(chDedupFiles,chSNPS_ref,ch_fasta,chSampleInfo,chIndexFiles)

    // Processo de SNP Fingerprint
    chSnpFingerprintComplete = snp_fingerprint(chDedupFiles, chSNPS_ref, ch_fasta, chSampleInfo, chIndexFiles).collect()

    // Processo SNP Footprint Clustering (executa apenas após a conclusão de snp_fingerprint para todas as amostras)
    snp_footprint_clustering(chSampleInfo,chRSNPFootprint,chSnpFingerprintComplete)

    enrichment(chEnrichmentScript,chDedupFiles,chSampleInfo)
    chFragDis = lenght_fragment_dist_step1(chDedupFiles,chSampleInfo)
    lenght_fragment_dist_step2(chRfrag_plotFragDist,chFragDis,chSampleInfo)

    chBWFiles = bedGraphToBigWig(chChromSizes,chPeakFiles,chSampleInfo)
    pileups_report(chSampleInfo,chChromSizes,chBWFiles,chPileUpBED,chRPileups)

    //Collect all files output and the pass to me program that will merge then
    //chAllFiles = chBWFiles.collectFile()
    //pileups_report_comp(chSampleDirPileUps,chChromSizes,chAllFiles,chPileUpBED,chRComparison)*/
}

