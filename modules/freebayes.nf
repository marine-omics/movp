process freebayes {

    publishDir "$params.outdir/freebayes"

    input:
    path(bam)
    path(bai)
    path fasta
    path fasta_fai
    path populations

    output:
    path("*.vcf.gz"), emit: vcf

    script:
    def prefix="freebayes"
    def chunksize= 10000
    def populations_file = populations.name != 'NO_FILE'  ? "--populations ${populations}" : ""

    def args = task.ext.args ?: ''

    if (task.cpus > 1) {
        """
        ls *.bam > bamlist.txt
        freebayes-parallel \\
            <(fasta_generate_regions.py $fasta_fai $chunksize) $task.cpus \\
            -f $fasta \\
            -L bamlist.txt \\
            $args \\
            $populations_file \\
            --strict-vcf > ${prefix}.vcf

        bgzip ${prefix}.vcf

        """

    } else {
        """
        ls *.bam > bamlist.txt
        freebayes \\
            -f $fasta \\
            -L bamlist.txt \\
            $args \\
            $populations_file \\
            --strict-vcf > ${prefix}.vcf

        bgzip ${prefix}.vcf

        """
    }

}


//--populations FILE
//                   Each line of FILE should list a sample and a population which
//                   it is part of.  The population-based bayesian inference model
//                   will then be partitioned on the basis of the populations.