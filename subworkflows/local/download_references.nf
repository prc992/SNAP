nextflow.enable.dsl=2

// Import the required processes from the modules
include { downloadGenome } from '../../modules/local/download'
include { createGenomeIndex } from '../../modules/local/createGenomeIndex'

workflow download_references {
    take:
        tuple val(genome), val(faGZFile), val(geneAnnotation), val(dacList), val(snp)
        path refDir

    //output:
    //    path "${genome}.fa", emit: genome
    //    path "${genome}.fa.*", emit: genome_index

    chGenome = downloadGenome(genome, faGZFile, geneAnnotation, dacList, snp, refDir)
    chGenomeIndex = createGenomeIndex(genome, faGZFile, geneAnnotation, dacList, snp, chGenome, refDir)

    emit: genome = chGenome
    emit: genome_index = chGenomeIndex

}