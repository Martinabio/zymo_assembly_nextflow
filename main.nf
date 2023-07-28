include {racon_polish_wf; map_reads_to_assembly} from './subworkflow/racon_polish.nf'

nextflow.preview.recursion=true

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

process sort_fastq {

    label 'process_high'
    label 'process_extreme_memory'

    container 'biocontainers/seqkit:2.5.0--h9ee0642_0'

    publishDir "${params.outdir}/"

    input:
    path in_fastq

    output:
    path "sorted_reads.fastq"

    script:
    """
    seqkit sort --reverse -l --threads $task.cpus ${in_fastq} > sorted_reads.fastq
    """

}

process assembly {

    label "process_high"
    label "process_extreme_memory"
    label "process_long"

    container 'biowilko/metamdbg:0.3'

    publishDir "${params.outdir}/unpolished_contigs/", mode: 'copy'

    input:
    path filtered_fastq
    output:
    path "metamdbg_out/contigs.fasta.gz"

    script:
    """
    /metaMDBG/build/bin/metaMDBG asm -t $task.cpus metamdbg_out/ ${filtered_fastq}
    """
}

process sort_convert_index_sam {

    label "process_medium"

    container 'biocontainers/samtools:1.3.1--h0cf4675_11'

    publishDir "${params.outdir}/sorted_bam/", mode: 'copy'

    input:
    path unsorted_sam

    output:
    path "sorted.bam"
    path "sorted.bam.bai"

    script:
    """
    samtools sort -O BAM -o sorted.bam --threads $task.cpus ${unsorted_sam} && samtools index sorted.bam
    """
}


process medaka_consensus {

    label "process_high"

    container 'biocontainers/medaka:1.8.0--py38hdaa7744_0'

    publishDir "${params.outdir}/polished_contigs/", mode: 'copy'

    input:
    path sorted_bam
    path sorted_bam_index
    path unpolished_contigs

    output:
    path "polished_contigs.fasta"

    script:
    """
    medaka consensus --threads $task.cpus --model ${params.medaka_model} ${sorted_bam} polished_contigs.hdf

    medaka stitch --threads $task.cpus polished_contigs.hdf ${unpolished_contigs} polished_contigs.fasta
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
        .collect()
        .set {combined_fastq_ch}

    // sort_fastq(combined_fastq_ch)

    assembly(combined_fastq_ch)

    assembly.out
        .collect()
        .set { assembly_ch}

    racon_polish_wf
        .recurse(assembly_ch, combined_fastq_ch)
        .times(4)

    map_reads_to_assembly(racon_polish_wf.out)

    sort_convert_index_sam(map_reads_to_assembly.out)

    medaka_consensus(sort_convert_index_sam.out, racon_polish_wf.out.contigs)
}

