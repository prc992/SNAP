process downloadSNPRef {
    label 'low_cpu_low_mem'
    tag "Dowloading - $genome" 
    publishDir "$path_sample_multiqc", mode : 'copy'

    container = params.containers.wget
    
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

    container = params.containers.wget
    
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

    container = params.containers.wget
    
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

process downloadGenome {

    label 'low_cpu_low_mem'
    tag "Dowloading - $genome" 
    //publishDir "$genomeOut", mode : 'copy'

    container = params.containers.wget
    
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