process mpileup_call {
//    publishDir "$params.outdir/bcftools", mode: 'copy'

    input:
    path(bam)
    path(bai)
    path fasta
    path fasta_fai
    val region

    output:
    path("*.vcf"), emit: vcf
//    path("*.vcf.gz.tbi"), emit: vcfi

    script:
    def callargs = task.ext.callargs ?: ''
    def args = task.ext.args ?: ''

    """
    ls *.bam > bamlist.txt
    bcftools mpileup --threads ${task.cpus} -r $region -Ou $args -f $fasta -b bamlist.txt | \\
    bcftools call --threads ${task.cpus} $callargs -v -Ov -o bcftools.${region}.vcf
    """

}

process mpileup_collect {

    publishDir "$params.outdir/bcftools", mode: 'copy'

    input:
    path(vcf)
    path(regions_file)

    output:
    path("*.vcf.gz"), emit: vcf
    path("*.vcf.gz.tbi"), emit: vcfi


    shell:
    '''
    while read region;do echo "bcftools.$region.vcf";done < !{regions_file} | xargs cat | vcffirstheader | vcfstreamsort -w 1000 > bcftools.vcf

    bgzip bcftools.vcf
    tabix bcftools.vcf.gz
    '''
}


process fasta_generate_chrs {

    input:
    path fasta
    path fasta_fai

    output:
    path("genome.chrs")

    script:

    def args = task.ext.args ?: ''

    """
    grep '>' ${fasta} | awk '{print \$1}' | sed 's/>//' > genome.chrs
    """

}

process gatk_gathervcfs {

    publishDir "$params.outdir/gatk", mode: 'copy'

    input:
      path(vcfs)

    output:
    path("gatk.vcf.gz"), emit: vcfz
    path("gatk.vcf.gz.tbi"), emit: vcfi    

    script:
    def args = task.ext.args ?: ''


    """
    bcftools concat ${vcfs} -O b -o gatk.vcf.gz
    tabix gatk.vcf.gz
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