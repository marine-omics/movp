
process faidx {

    input:
    path fasta

    output:
    path ("*.fai"), emit: fai

    script:
    """
    samtools faidx $fasta
    """
}

process sidx {

    input:
    path bam

    output:
    path("*.bai")

    script:
    """
    samtools index ${bam}
    """
}

// process samtools_idx {

//     input:
//     path bam

//     output:
//     tuple bam, path "*.bai" , emit: ibam

//     script:
//     """
//     samtools index $bam
//     """
// }
