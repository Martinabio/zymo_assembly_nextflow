process fastq_split {
    
    label 'process_single'

    container 'biocontainers/seqtk:1.4--he4a0461_1'

    input:
        path in_fastq

    output:
        path "*.fa"

    script:

    """
    seqtk split -n ${params.splits} split ${in_fastq}
    """

}

process fastp {
    
    label 'process_high'

    container 'biocontainers/fastp:0.23.4--hadf994f_1'

    input:
        path in_fastq

    output:
        path "${fastq_name}.filtered.fastq"

    script:

    fastq_name = in_fastq.getBaseName()

    """
    mv ${in_fastq} ${fastq_name}.fastq
    fastp --in1 ${fastq_name}.fastq --out1 ${fastq_name}.filtered.fastq --thread $task.cpus --low_complexity_filter -e ${params.min_mean_qual} -l ${params.min_length} 2> fastp.log
    """

}

process porechop {

    label 'process_high'

    container 'biowilko/porechop:0.1'

    input:
    path in_fastq
    output:
    path "trimmed.fastq"

    script:
    """
    porechop -i ${in_fastq} -o trimmed.fastq -t $task.cpus
    """
}

process flye_assembly {

    label "process_high"
    label "process_high_memory"
    label "process_long"

    container 'biocontainers/flye:2.9.2--py310h2b6aa90_2'

    publishDir "${params.outdir}/unpolished_contigs/", mode: 'copy'

    input:
    path filtered_fastq
    output:
    path "flye_output/assembly.fasta"

    script:
    """
    flye --nano-raw ${filtered_fastq} -o flye_output -g ${params.flye_genome_size} -t $task.cpus -i ${params.flye_iterations}
    """
}

process map_reads_to_assembly {

    label "process_medium"

    container 'biocontainers/minimap2:2.26--he4a0461_1'

    input:
    path flye_assembly
    path filtered_fastq

    output:
    path "unsorted.sam"

    script:
    """
    minimap2 -a -x map-ont -t $task.cpus ${flye_assembly} ${filtered_fastq} > unsorted.sam
    """
}

process sort_bam {

    label "process_medium"

    container 'biocontainers/samtools:1.3.1--h0cf4675_11'

    publishDir "${params.outdir}/sorted_bam/", mode: 'copy'

    input:
    path unsorted_sam

    output:
    path "sorted.bam"

    script:
    """
    samtools sort -O BAM -o sorted.bam --threads $task.cpus ${unsorted_sam}
    """
}

process medaka_consensus {

    label "process_high"

    container 'biocontainers/medaka:1.8.0--py38hdaa7744_0'

    publishDir "${params.outdir}/polished_contigs/", mode: 'copy'

    input:
    path sorted_bam

    output:
    path "polished_contigs.fasta"

    script:
    """
    medaka consensus --threads $task.cpus --model ${params.medaka_model} ${sorted_bam} polished_contigs.fasta
    """
}

workflow {
    Channel.of(file(params.input_fastq, type: "file", checkIfExists: true))
        .set {in_fastq_ch}

    if (params.split) {
        fastq_split(in_fastq_ch)
        fastq_split.out
            .flatten()
            .set { split_fastq_ch }
    } else {
        in_fastq_ch
            .set { split_fastq_ch }
    }

    fastp(split_fastq_ch)

    porechop(fastp.out)

    porechop.out
        .collectFile(name: "${params.outdir}/combined_trimmed_and_filtered.fastq", newLine: false)
        .set {combined_fastq_ch}

    flye_assembly(combined_fastq_ch)

    map_reads_to_assembly(flye_assembly.out, combined_fastq_ch)

    sort_bam(map_reads_to_assembly.out)

    medaka_consensus(sort_bam.out)
}

