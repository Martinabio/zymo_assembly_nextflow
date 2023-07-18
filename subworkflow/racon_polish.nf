include {map_reads_to_assembly} from '../main.nf'

process racon_polish {

    label "process_high"

    container 'biocontainers/racon:1.5.0--h21ec9f0_2'

    publishDir "${params.outdir}/racon_polished_contigs", mode: 'copy'

    input:
    path fastq_reads
    path sam_file
    path unpolished_contigs

    output:
    path "${fastq_reads}", emit: reads
    path "racon_polished_contigs.fasta", emit: polished_contigs

    script:
    """
    racon consensus --threads $task.cpus -u ${fastq_reads} ${sam_file} ${unpolished_contigs} > racon_polished_contigs.fasta 
    """
}

workflow racon_polish_wf {
    take: 
        contigs
        reads
    main:
        map_reads_to_assembly(contigs, reads)

        racon_polish(reads, map_reads_to_assembly.out, contigs)
    
    emit:
        contigs = racon_polish.out.polished_contigs
        reads = racon_polish.out.reads
}