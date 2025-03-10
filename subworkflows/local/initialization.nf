nextflow.enable.dsl=2

// Import the required processes from the modules
include {createSamplesheetFasta} from '../../modules/local/createSamplesheet'
include {fastqc} from '../../modules/local/fastqc'
include {multiqc} from '../../modules/local/multiqc'
include {moveSoftFiles} from '../../modules/local/moveSoftFiles'

workflow INITIALIZATION {


    take:
    chMultiQCConfig

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
    if (params.samplesheetfasta) {
        //println "Using provided samplesheet: ${params.samplesheet}"
        chSampleSheetFasta = Channel.fromPath(params.samplesheetfasta)
    } else {
        //println "Creating samplesheet because none was provided."
        chSampleSheetFasta = createSamplesheetFasta(
            params.sample_dir_fasta, 
            params.enrichment_mark ?: 'no_enrichment_mark'
        )
    }

    // Read the SampleSheet provided by the user or created by the pipeline
    chSampleInfo = chSampleSheetFasta \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.enrichment_mark, row.read1, row.read2) }
    
    // Run FastQC on the samples
    chFastaQC = fastqc(chSampleInfo)

    chFastaQCAll = chFastaQC.collect()

    // Filter only the files that will be used in the MultiQC report and remove duplicates
    chOnlyFiles = chFastaQCAll
        .map { values -> 
            values.findAll { 
                it instanceof Path && ( 
                    it.toString().endsWith(".yml") || 
                    it.toString().endsWith(".zip") || 
                    it.toString().endsWith(".html") )
            }
        }
        .flatten() // Garante que os arquivos estejam em um único fluxo
        .reduce( [:] as LinkedHashMap ) { acc, file -> 
            acc.putIfAbsent(file.getName(), file) // Mantém apenas a primeira ocorrência do nome do arquivo
            acc
        }
        .map { it.values().toList() } // 🔹 Converte para uma lista
        chFilesReportInitialization = chOnlyFiles.collect()

    chInitReport = multiqc(chFastaQCAll,chFilesReportInitialization,chMultiQCConfig)
    moveSoftFiles(chInitReport)

    emit: sample_info = chSampleInfo
    emit: genomes_info = chGenomesInfo
    emit: fastqc_files = chFastaQC
    emit: ref_dir = refDir
    emit: files_report_initialization = chFilesReportInitialization
    emit: init_report = chInitReport

}