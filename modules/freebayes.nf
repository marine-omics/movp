process freebayes_parallel {

    publishDir "$params.outdir/freebayes", mode: 'copy'

    input:
    path(bam)
    path(bai)
    path fasta
    path fasta_fai
    path populations

    output:
    path("*.vcf.gz"), emit: vcf
    path("*.vcf.gz.tbi"), emit: vcfi

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
        tabix ${prefix}.vcf.gz
        find ./ -type l | xargs rm
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
        tabix ${prefix}.vcf.gz
        find ./ -type l | xargs rm
        """
    }

}


process freebayes {

    input:
    path(bam)
    path(bai)
    path fasta
    path fasta_fai
    val region
    path populations

    output:
    path("*.vcf"), emit: vcf

    script:
    def prefix="freebayes"
    def populations_file = populations.name != 'NO_FILE'  ? "--populations ${populations}" : ""

    def args = task.ext.args ?: ''

    """
        ls *.bam > bamlist.txt
        freebayes --region ${region} \\
            -f $fasta \\
            -L bamlist.txt \\
            $args \\
            $populations_file \\
            --strict-vcf > ${prefix}.${region}.vcf
    """

}

process freebayes_collect {

    publishDir "$params.outdir/freebayes", mode: 'copy'

    input:
    path(vcf)
    path(regions_file)

    output:
    path("*.vcf.gz"), emit: vcf
    path("*.vcf.gz.tbi"), emit: vcfi


    shell:
    '''
    while read region;do echo "freebayes.$region.vcf";done < !{regions_file} | xargs cat | vcffirstheader | vcfstreamsort -w 1000 > freebayes.vcf

    bgzip freebayes.vcf
    tabix freebayes.vcf.gz
    '''
}


process fasta_generate_regions {

    input:
    path fasta
    path fasta_fai
    val chunksize

    output:
    path("*.regions")

    script:

    def args = task.ext.args ?: ''

    """
    fasta_generate_regions.py $fasta_fai $chunksize > freebayes.regions
    """


}



//--populations FILE
//                   Each line of FILE should list a sample and a population which
//                   it is part of.  The population-based bayesian inference model
//                   will then be partitioned on the basis of the populations.