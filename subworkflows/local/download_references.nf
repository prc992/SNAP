nextflow.enable.dsl=2

// Import the required processes from the modules
include {downloadGenome} from '../../modules/local/download'
include {downloadDACFile} from '../../modules/local/download'
include {downloadSNPRef} from '../../modules/local/download'
include {createGenomeIndex} from '../../modules/local/createGenomeIndex'
include {fetch_chrom_sizes} from '../../modules/local/fetch_chrom_sizes'

workflow DOWNLOAD_REFERENCES {

    take:
    chGenomesInfo
    chrefDir

    main:

    chGenome = downloadGenome(chGenomesInfo, chrefDir)
    chGenomeIndex = createGenomeIndex(chGenomesInfo,chGenome, chrefDir)
    chChromSizes = fetch_chrom_sizes(chGenomesInfo,chrefDir)
    chDACFileRef = downloadDACFile(chGenomesInfo,chrefDir)
    chSNPS_ref = downloadSNPRef(chGenomesInfo)

    emit: genome = chGenome
    emit: genome_index = chGenomeIndex
    emit: chrom_sizes = chChromSizes
    emit: dac_file_ref = chDACFileRef
    emit: snp_ref = chSNPS_ref

}