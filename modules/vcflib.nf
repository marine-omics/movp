process vcfwave {

    publishDir "$params.outdir/vcfwave", mode: 'copy'

    input:
    path(vcfin)
    path(vcfini)

    output:
    path("*.vcf.gz"), emit: vcf
    path("*.vcf.gz.tbi"), emit: vcfi

    script:
    """

    """
}