process signalIntensityCalculation {
    label 'low_cpu_low_mem'
    container = params.containers.r_data_analysis
    tag "Sample - $sampleId" 

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/signal_intensity/", mode: 'copy'

    input:
    tuple val(sampleId),val(enrichment_mark),val(control),val(read_method),path (bed_file),val (yml)
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
    # create an empty file to satisfy the output requirement if no enrichment mark is available
    touch ${sampleId}_${enrichment_mark}_Signal_Intensity_Matrix.csv
    """
}