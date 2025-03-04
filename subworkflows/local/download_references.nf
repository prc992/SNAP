nextflow.enable.dsl=2

// Import the required processes from the modules
include { downloadGenome } from '../../modules/local/download'
include { createGenomeIndex } from '../../modules/local/createGenomeIndex'

workflow DOWNLOAD_REFERENCES {

    take:
    chGenomesInfo
    chrefDir

    main:

    chGenome = downloadGenome(chGenomesInfo, chrefDir)
    chGenomeIndex = createGenomeIndex(chGenomesInfo,chGenome, chrefDir)

    emit: genome = chGenome
    emit: genome_index = chGenomeIndex

}