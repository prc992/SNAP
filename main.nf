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
    
    input:
    path qc_reports_dir

    exec:
    path_sample_multiqc =  params.output_dir + "/reports/multiqc/" 

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc ${qc_reports_dir} -o ./
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

    output:
    path "snap-samplesheet-*.csv"

    script:
    """
    now=\$(date +'%Y-%m-%d-%H-%M-%S')
    filename="snap-samplesheet-\$now.csv"
    echo "sampleId,path,read1,read2" > \$filename

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
        echo "\$sampleId,${output_dir},\$read1,\$read2" >> \$filename
    done
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
    println "SNAP: Streamlined Nextflow Analysis Pipeline for profiling circulating histone modifications identifies tumor epigenomic signatures in cancer plasma"
    println "GitHub repository: ${githubPath}"
    println "Release version: ${releaseVersion}"


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
    chGeneAnotation = downloadGeneAnotation(params.genome,refDir)
    chChromSizes = fetch_chrom_sizes(params.genome,refDir)
    chDACFileRef = downloadDACFile(params.genome,refDir)
    
    // Create the output directory if it doesn't exist
    """
    mkdir -p ${projectDir}/${params.output_dir}
    """.execute().waitFor()

    chSampleSheet = createSamplesheet(params.sample_dir, params.output_dir)

    chSampleInfo = chSampleSheet \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,"${projectDir}/${row.path}", row.read1, row.read2) }

    fastqc(chSampleInfo)
    chTrimFiles = trim(chSampleInfo)
    chAlignFiles = align(chTrimFiles,chGenome,chGenomeIndex)    
    chSortedFiles = sort_bam(chAlignFiles)
    lib_complex(chSortedFiles)
    chUniqueFiles = unique_sam(chSortedFiles)
    chDedupFiles = dedup(chUniqueFiles)
    chDACFilteredFiles = dac_exclusion(chDedupFiles,chDACFileRef)
    chIndexFiles = index_sam(chDACFilteredFiles)

    chPeakFiles = peak_bed_graph(chDACFilteredFiles)
    uropa(chPeakFiles,chGeneAnotation)
    chBedFiles = bam_to_bed(chDACFilteredFiles)
    unique_frags(chBedFiles)


    // Processo de SNP Fingerprint
    chSnpFingerprintComplete = snp_fingerprint(chIndexFiles, chSNPS_ref, chGenome).collect()
    //chSnpFingerprintComplete = snp_fingerprint(chIndexFiles, chSNPS_ref, chGenome,chGenomeIndex).collect()

    // Ver Depois
    // Processo SNP Footprint Clustering (executa apenas após a conclusão de snp_fingerprint para todas as amostras)
    //snp_footprint_clustering(chSnpFingerprintComplete,chRSNPFootprint)

    enrichment(chDACFilteredFiles,chEnrichmentScript)
    chFragDis = lenght_fragment_dist_step1(chDACFilteredFiles)
    lenght_fragment_dist_step2(chFragDis,chRfrag_plotFragDist)

    chBWFiles = bedGraphToBigWig(chPeakFiles,chChromSizes)
    pileups_report(chBWFiles,chChromSizes,chPileUpBED,chRPileups)

    // Collect QC reports from FastQC, Trimmed files, and other relevant steps
    chAllQCFiles = Channel.from(
        fastqc.out.collectFile(),
        trim.out.collectFile(),
        align.out.collectFile()
        // Add other relevant QC output channels if applicable
    )

    // Pass the QC files to the MultiQC process
    multiqc(chAllQCFiles)


    /*//Collect all files output and the pass to me program that will merge then
    //chAllFiles = chBWFiles.collectFile()
    //pileups_report_comp(chSampleDirPileUps,chChromSizes,chAllFiles,chPileUpBED,chRComparison)*/
}

