process fragLenHist {
    label 'med_cpu_med_mem'
    container = params.containers.python
    tag "All Samples"

    publishDir "$path_sample_frags", mode : 'copy'

    input:
    path (raw_fragments)
    each path (frag_len_header_multiqc)
    each path (chCalcFragHist)
    tuple val(sampleId), val(enrichment_mark),val(path_analysis),val(read1), val(read2)

    exec:
    path_sample_frags = path_analysis + "/frag/"

    output:
    tuple path ("frag_len_hist.txt"),path ("frag_len_mqc.yml")

    script:
    """
    python $chCalcFragHist \\
        --frag_path "*fragment_sizes.txt" \\
        --output frag_len_hist.txt

    cat $frag_len_header_multiqc frag_len_hist.txt > frag_len_mqc.yml
    """

}