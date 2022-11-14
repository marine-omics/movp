

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
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bai")

    script:
    """
    samtools index $bam
    """
}


process flagstat {

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path ("*.flagstat"), emit: flagstat

    script:

    def outfile = "${meta.sample}.flagstat"

    """
    samtools flagstat $bam > $outfile
    """
}

process stat {

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path ("*.stat"), emit: stat

    script:

    def outfile = "${meta.sample}.stat"

    """
    samtools stat $bam > $outfile
    """
}

process idxstat {

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path ("*.idxstat"), emit: idxstat

    script:

    def outfile = "${meta.sample}.idxstat"

    """
    samtools idxstats $bam > $outfile
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
