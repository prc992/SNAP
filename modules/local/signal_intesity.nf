process signalIntensityCalculation {
    label 'low_cpu_low_mem'
    container = params.containers.r_signal_intensity
    tag "Sample - $sampleId" 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/signal_intensity/", mode: 'copy'

    input:
    tuple val(sampleId),val(enrichment_mark),val(control),val(read_method),path (bed_file),val (yml)
    each path(chDACFileRef)
    each path(chRMEDIPSignalCalculation)
    each path(chRMARKSSignalCalculation)
    each path(chRegions_of_interest_MEDIP_signal)
    each path(chRegions_of_interest_MARKS_signal)
    each path(chHousekeeping_MEDIP_signal)
    each path(chHousekeeping_H3K4ME3_signal)
    each path(chHousekeeping_H3K27AC_signal)

    output:
    path ("*_Signal_Intensity_Matrix.csv")

    script:
    """
    if [ "\$(echo '${enrichment_mark}' | tr '[:upper:]' '[:lower:]')" = "medip" ]; then
        Rscript ${chRMEDIPSignalCalculation} ${bed_file} ${chRegions_of_interest_MEDIP_signal} ${chHousekeeping_MEDIP_signal}
    elif [ "\$(echo '${enrichment_mark}' | tr '[:upper:]' '[:lower:]')" = "h3k4me3" ]; then
        Rscript ${chRMARKSSignalCalculation} ${bed_file} ${chDACFileRef} ${chRegions_of_interest_MARKS_signal} ${chHousekeeping_H3K4ME3_signal} ${enrichment_mark}
    elif [ "\$(echo '${enrichment_mark}' | tr '[:upper:]' '[:lower:]')" = "h3k27ac" ]; then
        Rscript ${chRMARKSSignalCalculation} ${bed_file} ${chDACFileRef} ${chRegions_of_interest_MARKS_signal} ${chHousekeeping_H3K27AC_signal} ${enrichment_mark}
    fi
    
    # Verificar se o arquivo de saída foi criado, caso contrário criar um arquivo vazio
    if [ ! -f ${sampleId}_${enrichment_mark}_Signal_Intensity_Matrix.csv ]; then
        touch ${sampleId}_${enrichment_mark}_Signal_Intensity_Matrix.csv
    fi
    """
}