process createSamplesheet {
    label 'low_cpu_low_mem'
    tag "Creating Samplesheet" 

    publishDir "${workflow.projectDir}/${params.outputFolder}", mode : 'copy'

    input:
    val sample_dir
    val enrichment_mark

    output:
    path "snap-samplesheet-*.csv"

    script:
    """
    now=\$(date +'%Y-%m-%d-%H-%M-%S')
    filename="snap-samplesheet-\$now.csv"
    echo "sampleId,enrichment_mark,read1,read2" > \$filename

    for subfolder in \$(find ${sample_dir} -mindepth 1 -maxdepth 1 -type d); do
        sampleId=\$(basename \$subfolder)
        # Exclude hidden directories
        if [[ \$sampleId == .* ]]; then
            continue
        fi
        files=(\$(find \$subfolder -type f \\( -name '*.fq.gz' -o -name '*.fq' -o -name '*.fastq.gz' \\) | sort))
        read1=\$(realpath \${files[0]})
        read2=""
        if [ \${#files[@]} -gt 1 ]; then
            read2=\$(realpath \${files[1]})
        fi
        echo "\$sampleId,${enrichment_mark},\$read1,\$read2" >> \$filename
    done
    """
}