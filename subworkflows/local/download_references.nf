nextflow.enable.dsl=2

// Import the required processes from the modules
include {downloadGenome} from '../../modules/local/download'
include {downloadDACFile} from '../../modules/local/download'
include {downloadSNPRef} from '../../modules/local/download'
include {createGenomeIndex} from '../../modules/local/createGenomeIndex'
include {fetch_chrom_sizes} from '../../modules/local/fetch_chrom_sizes'
include {createSamplesheet} from '../../modules/local/createSamplesheet'

workflow DOWNLOAD_REFERENCES {

    take:
    chGenomesInfo
    chrefDir

    main:

    chGenome = downloadGenome(chGenomesInfo, chrefDir)
    chGenomeIndex = createGenomeIndex(chGenomesInfo,chGenome, chrefDir)
    chChromSizes = fetch_chrom_sizes(chGenomesInfo,chrefDir)
    chDACFileRef = downloadDACFile(chGenomesInfo,chrefDir)

     // If the 'samplesheet' parameter is provided, use it directly; otherwise, create a new samplesheet
    if (params.samplesheet) {
        //println "Using provided samplesheet: ${params.samplesheet}"
        chSampleSheet = Channel.fromPath(params.samplesheet)
    } else {
        //println "Creating samplesheet because none was provided."
        chSampleSheet = createSamplesheet(
            params.sample_dir, 
            params.enrichment_mark ?: 'no_enrichment_mark'
        )
    }

    // Read the SampleSheet provided by the user or created by the pipeline
    chSampleInfo = chSampleSheet \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.enrichment_mark, row.read1, row.read2) }

    chSNPS_ref = downloadSNPRef(chGenomesInfo)

    emit: genome = chGenome
    emit: genome_index = chGenomeIndex
    emit: chrom_sizes = chChromSizes
    emit: dac_file_ref = chDACFileRef
    emit: sample_info = chSampleInfo
    emit: snp_ref = chSNPS_ref

}