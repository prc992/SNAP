process createSMaSHFingerPrint{
    label 'med_cpu_med_mem'

    container = params.containers.python

    tag "All Samples"   

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/SMaSH/", mode : 'copy'
    
    input:
    path (chSNPSMaSH)
    path (chSNPS_ref)
    path (chBamAndBai)

    output:
    path('pval_out.txt')

    script:
    """
    python3 $chSNPSMaSH -i $chSNPS_ref ALL
    """
}

process createSMaSHFingerPrintPlot{
    label 'med_cpu_med_mem'

    container = params.containers.snp_smash_plot

    tag "All Samples"   

    publishDir "${workflow.projectDir}/${params.outputFolder}/reports/SMaSH/", mode : 'copy'
    
    input:
    path (chSMaSHOutout)
    path (chSNPSMaSHPyPlot)

    output:
    path('*.jpg')

    script:
    """
    python3 $chSNPSMaSHPyPlot -i $chSMaSHOutout -o Dendrogram_of_Samples_by_SNP_Profile.jpg
    """
}