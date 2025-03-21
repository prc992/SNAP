process createSamplesheetFasta {
    label 'low_cpu_low_mem'
    tag "Creating Samplesheet Fasta" 

    publishDir "${workflow.projectDir}/${params.outputFolder}", mode : 'move'

    input:
    val sample_dir
    val enrichment_mark
    val sample_control

    output:
    path "snap-samplesheet-fasta-*.csv"

    script:
    """
    now=\$(date +'%Y-%m-%d-%H-%M-%S')
    filename="snap-samplesheet-fasta-\$now.csv"
    echo "sampleId,enrichment_mark,read1,read2,control" > \$filename

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

        if [ "\$sampleId" == "${sample_control}" ]; then
            echo "\$sampleId,${enrichment_mark},\$read1,\$read2," >> \$filename
        else
            echo "\$sampleId,${enrichment_mark},\$read1,\$read2,${sample_control}" >> \$filename
        fi
        
    done
    """
}

process createSamplesheetBam {
    label 'low_cpu_low_mem'
    tag "Creating Samplesheet BAM" 

    publishDir "${workflow.projectDir}/${params.outputFolder}", mode : 'copy'

    input:
    val sample_dir
    val enrichment_mark

    output:
    path "snap-samplesheet-bam-*.csv"

    script:
    """
    now=\$(date +'%Y-%m-%d-%H-%M-%S')
    filename="snap-samplesheet-bam-\$now.csv"
    echo "sampleId,enrichment_mark,bam,control" > \$filename

    for subfolder in \$(find ${sample_dir} -mindepth 1 -maxdepth 1 -type d); do
        sampleId=\$(basename \$subfolder)
        # Exclude hidden directories
        if [[ \$sampleId == .* ]]; then
            continue
        fi
        files=(\$(find \$subfolder -type f \\( -name '*.bam' \\) | sort))
        bam=\$(realpath \${files[0]})


        if [ "\$sampleId" == "${sample_control}" ]; then
            echo "\$sampleId,${enrichment_mark},\$bam" >> \$filename
        else
            echo "\$sampleId,${enrichment_mark},\$bam$,${sample_control}" >> \$filename
        fi
            
    done
    """
}