nextflow.enable.dsl=2

// Import the required processes from the modules
include {createSamplesheet} from '../../modules/local/createSamplesheet'
include {fastqc} from '../../modules/local/fastqc'
include {multiqc_initialization} from '../../modules/local/multiqc'

workflow INITIALIZATION {

    main:

    // Create the genome directory if it doesn't exist
    """
    mkdir -p ${projectDir}/ref_files/genome
    """.execute().waitFor()
    refDir = Channel.fromPath("${projectDir}/ref_files/genome")

    // Read the GenomePaths spreadsheet and filter the row that matches the genome parameter
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
    
    // Run FastQC on the samples
    chFastaQC = fastqc(chSampleInfo)

    chFastaQCAll = chFastaQC.collect()

    multiqc_initialization(chFastaQCAll)

    emit: sample_info = chSampleInfo
    emit: genomes_info = chGenomesInfo
    emit: fastqc_files = chFastaQC
    emit: ref_dir = refDir

}