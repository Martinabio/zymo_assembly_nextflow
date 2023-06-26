process porechop {

    label 'process_medium'
    label 'process_high_memory'

    container 'biowilko/zymo_assembly_nextflow:0.0.1'

    publishDir "${params.outdir}/trimmed_reads/", mode: 'copy'

    input:
    path in_fastq
    output:
    path "trimmed.fastq"

    script:
    """
    porechop -i ${in_fastq} -o trimmed.fastq -t $task.cpus
    """
}

process filtlong {

    label "process_low"

    container 'biowilko/zymo_assembly_nextflow:0.0.1'

    publishDir "${params.outdir}/filtered_reads/", mode: 'copy'

    input:
    path trimmed_fastq
    output:
    path "filtered.fastq"

    """
    filtlong --min-length ${params.min_length} --min_mean_q ${params.min_mean_qual} ${trimmed_fastq} > filtered.fastq
    """

}

process flye_assembly {

    label "process_high"

    container 'biowilko/zymo_assembly_nextflow:0.0.1'

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

    container 'biowilko/zymo_assembly_nextflow:0.0.1'

    input:
    path flye_assembly
    path filtered_fastq

    output:
    path "sorted.bam"

    script:
    """
    minimap2 -a -x map-ont -t $task.cpus ${flye_assembly} ${filtered_fastq} | samtools sort -o sorted.bam
    """
}

process medaka_consensus {

    label "process_high"

    container 'biowilko/zymo_assembly_nextflow:0.0.1'

    publishDir "${params.outdir}/polished_contigs/", mode: 'copy'

    input:
    path sorted_bam

    output:
    path "polished_contigs.fasta"

    script:
    """
    medaka polish --threads $task.cpus --model ${params.medaka_model} ${sorted_bam} polished_contigs.fasta
    """
}

workflow {
    Channel.of(file(params.input_fastq, type: "file", checkIfExists: true))
        .set {in_fastq_ch}
    
    porechop(in_fastq_ch)

    filtlong(porechop.out)

    flye_assembly(filtlong.out)

    map_reads_to_assembly(flye_assembly.out, filtlong.out)

    medaka_consensus(map_reads_to_assembly.out)
}

