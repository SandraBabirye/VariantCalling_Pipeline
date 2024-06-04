nextflow.enable.dsl=2

process bwa_index {
    publishDir params.outdir, mode: 'copy'

    input:
    path ref

    output:
    path("${ref}.bwt")

    script:
    """
    if [ ! -f "${ref}.bwt" ]; then
        bwa index $ref
    fi
    """
}

process bwa_mem {
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_R1)
    path ref
    path indices

    output:
    path 'SEaligned.sam'

    script:
    """
    bwa mem -t 16 -aM $ref -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" $fastq_R1 > SEaligned.sam
    """
}



