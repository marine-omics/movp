process gatk4_createsequencedict {

    input:
    path fasta

    output:
    path "*.dict" , emit: dict

    script:
    """
    gatk CreateSequenceDictionary -REFERENCE $fasta -OUTPUT ${fasta.baseName}.dict
    """
}

process gatk4_createintervallist {
  input:
    path(reffaidx)
    path(refdict)

  output:
    path("*.interval_list")


  script:
    def outfile = "${refdict.baseName}.interval_list"

  """
  cat $reffaidx | awk '{OFS="\t";print \$1,0,\$2}' > ref.bed
  gatk BedToIntervalList -I ref.bed -O $outfile -SD $refdict
  """
}

process gatk_scatterintervals {
  input:
    path(intervallist)
    val(chunksize)

  output:
    path("scatter_temp*.interval_list")

  script:

  """
  mkdir scatter
  gatk IntervalListTools -I $intervallist --SCATTER_CONTENT $chunksize -O scatter
  for f in \$(find scatter -name '*.interval_list');do 
    tn=\$(echo \$f | sed 's:/:_:g');
    cp \$f \$tn
  done
  """
}


process fastq2ubam {

    input:
      tuple val(meta), path(reads)

    output:
      tuple val(meta), path("*.bam"), emit: ubam
//      val(meta), emit: meta
  //    path "*.bam", emit: bam

    script:

      def read_group  = "${meta.sample}.${meta.flowcell}.${meta.lane}"
      if (meta.single_end) {
        """
        gatk --java-options "-Xmx${task.memory.giga}G" FastqToSam \\
        -FASTQ ${reads[0]} \\
          -OUTPUT ${read_group}.bam \\
          -READ_GROUP_NAME ${read_group} \\
          -SAMPLE_NAME ${meta.sample} \\
          -LIBRARY_NAME ${meta.sample} \\
          -PLATFORM_UNIT ${read_group} \\
          -PLATFORM ILLUMINA
        """
        } else {
        """
        gatk --java-options "-Xmx${task.memory.giga}G" FastqToSam \\
        -FASTQ ${reads[0]} \\
        -FASTQ2 ${reads[1]} \\
          -OUTPUT ${read_group}.bam \\
          -READ_GROUP_NAME ${read_group} \\
          -SAMPLE_NAME ${meta.sample} \\
          -LIBRARY_NAME ${meta.sample} \\
          -PLATFORM_UNIT ${read_group} \\
          -PLATFORM ILLUMINA
        """
        }
}


process markadapters {

    input:
      tuple val(meta), path(ubam)
//      val(meta)
  //    path(ubam)

    output:
      tuple val(meta), path("*marked.bam"), emit: mdbam
//      val(meta), emit: meta
  //    path "*marked.bam", emit: mdbam

    script:

    def args = task.ext.args ?: ''

    def read_group  = "${meta.sample}.${meta.flowcell}.${meta.lane}"

    """
    gatk --java-options "-Xmx${task.memory.giga}G" MarkIlluminaAdapters \
    -I $ubam \
    -M ${read_group}_txt \
    $args \
    -O ${read_group}_marked.bam 
    """
}



process bwa_mem_gatk {

    input:
    tuple val(meta), path(mkbam), path(ubam)
    path(genome)
    path(index)
    path(dict)

    output:
    tuple val(meta), path("*mapped.bam"), emit: bam

    script:

    def halfmem = task.memory.giga/2
    def threads = task.cpus*2

    def outfile = "${mkbam.baseName}_mapped.bam"

    """
    gatk --java-options "-Xmx${halfmem}G" SamToFastq \
    -I $mkbam \
    -FASTQ /dev/stdout \
    -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
    |  \
    bwa mem -M -t ${threads} -p $genome /dev/stdin \
    | \
    gatk --java-options "-Xmx${halfmem}G" MergeBamAlignment \
    -ALIGNED_BAM /dev/stdin \
    -UNMAPPED_BAM $ubam \
    -OUTPUT $outfile \
    -R $genome -CREATE_INDEX true -ADD_MATE_CIGAR true \
    -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true \
    -INCLUDE_SECONDARY_ALIGNMENTS true -MAX_INSERTIONS_OR_DELETIONS -1 \
    -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS 
    """
}


process gatk_mark_duplicates {

  publishDir "$params.outdir/mapped_marked_bams", mode: 'copy'  

  input:
    tuple val(sample), path(bam)

  output:
    tuple val(sample), path("*.bam"), emit: mbam
    tuple val(sample), path("*_duplicatemetrics.txt"), emit: dupmetrics

  script:

  def args = task.ext.args ?: ''

  def outfile = "${bam.baseName}_marked.bam"

  """
    gatk --java-options "-Xmx${task.memory.giga}G" MarkDuplicates \
      -I $bam -O $outfile \
      $args \
      -METRICS_FILE ${bam.baseName}_duplicatemetrics.txt
  """
}


process gatk_haplotype_caller {

    input:
      tuple val(sample), path(bam), path(bai), path(interval)
      path(genome)
      path(index)
      path(dict)


    output:
    tuple val(sample),val(interval.baseName), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), emit: gvcf

    script:
    def args = task.ext.args ?: ''
    def outfile = "${bam.baseName}.${interval.baseName}.g.vcf.gz"
    def threads = task.cpus*2

    """

    gatk --java-options "-Xmx${task.memory.giga}G" \
        HaplotypeCaller \
        --native-pair-hmm-threads ${threads} \
        -R $genome \
        -I $bam \
        -L $interval \
        -O $outfile \
        -contamination 0 -ERC GVCF
    """

}


process gatk_genomicsdb_import {

    input:
      tuple val(regionid), path(intervallist), path(samplemap)

    output:
    tuple val(regionid), path(intervallist), path("*.gatk.db"), emit: gdb

    script:
    def args = task.ext.args ?: ''

    """
    gatk --java-options "-Xmx${task.memory.giga-1}g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ${regionid}.gatk.db \
        --sample-name-map $samplemap \
        --reader-threads ${task.cpus} \
        --batch-size 50 \
        --L ${intervallist}
    chmod a+rx *.db
    """

}




process gatk_genotypegvcfs {

    input:
      tuple val(regionid), path(intervallist), path(gdb)
      path(genome)
      path(index)
      path(dict)

    output:
      tuple val(regionid), path("${regionid}.vcf.gz"), emit: vcfz

    script:
    def args = task.ext.args ?: ''
    def heap_gb = (task.memory.giga * params.java_heap_fraction) as int

    """
    gatk --java-options "-Xmx${heap_gb}g -Xms${heap_gb}g" \
        GenotypeGVCFs \
        -R $genome \
        -O ${regionid}.vcf.gz \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        -V gendb://${gdb} \
        -L $intervallist
    """

}


process gatk_mergevcfs {

  publishDir "$params.outdir/gatk", mode: 'copy'

  input:
    path(vcflist)
    path(dict)

  output:
    path("gatk.vcf.gz"), emit: vcfz
    path("gatk.vcf.gz.tbi"), emit: vcfi 

  script:

    """
    gatk --java-options "-Xmx{task.memory.giga}g -Xms{task.memory.giga}g" \
        MergeVcfs \
        -I $vcflist \
        -D $dict \
        -O gatk.vcf.gz
    """
}
