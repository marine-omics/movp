process freebayes {

    publishDir "$params.outdir/freebayes"

    input:
    path(bam)
    path(bai)
    path fasta
    path fasta_fai

    output:
    path("*.vcf.gz"), emit: vcf

    script:
    def prefix="freebayes"
    def chunksize= 10000

    if (task.cpus > 1) {
        """
        ls *.bam > bamlist.txt
        freebayes-parallel \\
            <(fasta_generate_regions.py $fasta_fai $chunksize) $task.cpus \\
            -f $fasta \\
            -L bamlist.txt \\
            -E -1 -m 30 -q 20 -K --strict-vcf > ${prefix}.vcf

        bgzip ${prefix}.vcf

        """

    } else {
        """
        ls *.bam > bamlist.txt
        freebayes \\
            -f $fasta \\
            -L bamlist.txt \\
            -E -1 -m 30 -q 20 -K --strict-vcf > ${prefix}.vcf

        bgzip ${prefix}.vcf

        """
    }

}
