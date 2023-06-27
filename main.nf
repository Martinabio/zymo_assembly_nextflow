process fastp {
    
    label 'process_medium'

    publishDir "${params.outdir}/fastp/", mode: 'copy'

    container 'biowilko/zymo_assembly_nextflow:0.0.3'

    input:
        path in_fastq

    output:
        path "trimmed.fastq.gz"

    script:

    """
    fastp --in1 ${fastq} --out1 trimmed.fastq.gz --thread $task.cpus --low_complexity_filter -e ${params.min_mean_qual} -l ${params.min_length} 2> fastp.log
    """

}

process flye_assembly {

    label "process_high"

    container 'biowilko/zymo_assembly_nextflow:0.0.3'

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

    container 'biowilko/zymo_assembly_nextflow:0.0.3'

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

    container 'biowilko/zymo_assembly_nextflow:0.0.3'

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
    
    pychopper(in_fastq_ch)

    flye_assembly(pychopper.out)

    map_reads_to_assembly(flye_assembly.out, pychopper.out)

    medaka_consensus(map_reads_to_assembly.out)
}

