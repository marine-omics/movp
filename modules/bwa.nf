process bwa_index {

    input:
    path fasta

    output:
    path "*" , emit: index

    script:
    """
    bwa index $fasta
    """
}
