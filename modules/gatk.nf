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
        gatk FastqToSam \\
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
        gatk FastqToSam \\
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

    def read_group  = "${meta.sample}.${meta.flowcell}.${meta.lane}"

    """
    gatk MarkIlluminaAdapters \
    -I $ubam \
    -M ${read_group}_txt \
    -O ${read_group}_marked.bam 

    """
}



process bwa_mem_gatk {

    publishDir "$params.outdir/mapped_bams"

    input:
    tuple val(meta), path(mkbam), path(ubam)
    path(genome)
    path(index)
    path(dict)

    output:
    tuple val(meta), path("*mapped.bam"), emit: bam
//    path "*mapped.bam", emit: bam

    script:

    def outfile = "${mkbam.baseName}_mapped.bam"

    """
    gatk SamToFastq \
    -I $mkbam \
    -FASTQ /dev/stdout \
    -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \
    |  \
    bwa mem -M -t $task.cpus -p $genome /dev/stdin \
    | \
    gatk MergeBamAlignment \
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

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path("*.bam"), emit: mbam


  script:

  def outfile = "${bam.baseName}_marked.bam"

  """
    gatk MarkDuplicates \
      -I $bam -O $outfile \
      -METRICS_FILE ${bam.baseName}_markduplicates_txt
  """
}


