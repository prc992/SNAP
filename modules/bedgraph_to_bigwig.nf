process bedgraph_to_bigwig {
    // Define os arquivos de entrada e saída

    container = '4dndcic/4dn-bedgraphtobigwig:v6'

    publishDir "${workflow.projectDir}/${params.outputFolder}/peaks/${sampleId}", mode : 'copy'

    input:
    tuple val(sampleId),path(chbedgraph),val (_)
    each path (chrom_sizes)

    output:
    tuple path ("*.bw"),path ("bedgraph_to_bigwig_mqc_versions.yml")

    exec:
    str_bw = sampleId + '.bw'

    script:
    """
    bedGraphToBigWig $chbedgraph $chrom_sizes $str_bw

    cat <<-END_VERSIONS > bedgraph_to_bigwig_mqc_versions.yml
    "${task.process}":
        ucsc: $params.containers.bedgraphtobigwig_version
    END_VERSIONS
    """
}