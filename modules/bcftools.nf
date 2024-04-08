process mpileup_call {
    publishDir "$params.outdir/bcftools", mode: 'copy'

    input:
    path(bam)
    path(bai)
    path fasta
    path fasta_fai

    output:
    path("*.vcf.gz"), emit: vcf
    path("*.vcf.gz.tbi"), emit: vcfi

    script:
    def args = task.ext.args ?: ''

    """
    bcftools mpileup -Ou $args -f $fasta $bam | \\
    bcftools call --threads $task.cpus -mv -Oz -o bcftools.vcf.gz

    tabix bcftools.vcf.gz
    """

}



process gatk_gathervcfs {

    publishDir "$params.outdir/gatk", mode: 'copy'

    input:
      path(vcfs)

    output:
    path("gatk.vcf.gz"), emit: vcfz

    script:
    def args = task.ext.args ?: ''


    """
    bcftools concat ${vcfs} -O b -o gatk.vcf.gz
    """

}



// process isec_union {
//   publishDir "$params.outdir/vcfunion", mode: 'copy'

//   input:
//   path(vcf)
//   path(vcfi)
//   val(n)

//   output:
//   path('union.vcf.gz'), emit:vcf
//   path('union.vcf.gz.tbi'), emit: vcfi

//   script:
//   """
//   bcftools isec $vcf -n $n -Oz  -o union.vcf.gz
//   tabix union.vcf.gz
//   """
// }