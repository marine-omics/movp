nextflow.enable.dsl=2

include { fastqc } from './modules/fastqc.nf'
include { multiqc_fastqc; multiqc_fastp; multiqc_bams } from './modules/multiqc.nf'
include { fastp } from './modules/fastp.nf'
include { fastq2ubam; markadapters; bwa_mem_gatk; gatk4_createsequencedict; gatk4_createintervallist; gatk_scatterintervals; gatk_mark_duplicates; gatk_haplotype_caller; gatk_genomicsdb_import; gatk_genotypegvcfs; gatk_mergevcfs } from './modules/gatk.nf'
include { bwa_index } from './modules/bwa.nf'
include { sidx; faidx; flagstat; stat; idxstat; samtools_merge } from './modules/samtools.nf'
include { freebayes; fasta_generate_regions; freebayes_collect } from './modules/freebayes.nf'
include { mpileup_call; mpileup_collect; fasta_generate_chrs; gatk_gathervcfs } from './modules/bcftools.nf'
include { name_by_sample } from './modules/util.nf'

workflow qc {
  take:
    fastqin
  main:
    fastqin | fastqc | collect | multiqc_fastqc
}

workflow preprocess {
  take:
    fastqin

  main:
    fastp(fastqin) 
    fastp.out.json | collect | multiqc_fastp

  emit:
    fastp.out.reads
}


workflow gatk_map {
  take:
    preads
    genome_fasta
    genome_index
    genome_dict

  main:
    ch_ubams = preads | fastq2ubam
    ch_marked_bams =  ch_ubams | markadapters
    ch_merge_bams = ch_ubams.join(ch_marked_bams)

    mapped_bams = bwa_mem_gatk(ch_merge_bams,genome_fasta,genome_index, genome_dict)

    mapped_bams.map { m,b -> [m.sample,b] }.groupTuple().branch{
        //Here there actually is a list, so size() works
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{bam_to_merge} 


    ch_merged_multi_bams = samtools_merge(bam_to_merge.multiple)
    ch_rename_single_bams = name_by_sample(bam_to_merge.single)

    ch_persample_bams = ch_merged_multi_bams.mix(ch_rename_single_bams)

    mapped_marked_bams = gatk_mark_duplicates(ch_persample_bams)

  emit:
    mapped_marked_bams
}

workflow bam_qc {
  take:
    bams

  main:
    fs = bams | flagstat | collect
    s = bams | stat | collect
    idxf = bams | idxstat | collect
    multiqc_bams(fs,s,idxf)
}



workflow call_variants {
  take:
    callers
    indexed_bams
    genome_fasta
    genome_fai
    genome_dict

  main:

    ch_bamcollection = indexed_bams.map{m,b,i -> b} | collect
    ch_baicollection = indexed_bams.map{m,b,i -> i} | collect

    if ( callers.contains('freebayes') ){
        regions_file = fasta_generate_regions(genome_fasta,genome_fai,params.fb_chunksize)
        ch_regions = regions_file.splitText().map{it -> it.trim()}


      // Freebayes
      ch_chunk_vcfs = freebayes(ch_bamcollection,ch_baicollection,genome_fasta,genome_fai,ch_regions,file(params.populations)) | collect

      freebayes_collect(ch_chunk_vcfs,regions_file)

    } 

    if ( callers.contains('bcftools') ){
      chrs_file = fasta_generate_chrs(genome_fasta,genome_fai)
      ch_chrs = chrs_file.splitText().map{it -> it.trim()}

      // bcftools
      ch_chr_vcfs = mpileup_call(ch_bamcollection,ch_baicollection,genome_fasta,genome_fai,ch_chrs) | collect  
      mpileup_collect(ch_chr_vcfs,chrs_file)
    } 

    if ( callers.contains('gatk') ){
      genome_intervals = gatk4_createintervallist(genome_fai,genome_dict)
      ch_gatk_scatter_intervals = gatk_scatterintervals(genome_intervals,params.gatk_chunksize) | flatten

      ch_hc_inputs = indexed_bams.combine(ch_gatk_scatter_intervals)


      // gatk
      ch_gvcfs = gatk_haplotype_caller(ch_hc_inputs,genome_fasta,genome_fai,genome_dict)

      samples_gvcfs_file = ch_gvcfs.collectFile{ 
        s,i,b,tbi -> [ "${i}.sample_map" , "$s\t$b\n" ]
      }.map {
        it -> [file(it).baseName, it ]
      }

      ch_sample_maps_per_region = ch_gatk_scatter_intervals.map{ 
        it -> [file(it).baseName, it ]
      }.join(samples_gvcfs_file)


      ch_gdb = gatk_genomicsdb_import(ch_sample_maps_per_region)

      ch_gatk_vcfs = gatk_genotypegvcfs(ch_gdb,genome_fasta,genome_fai,genome_dict) 

      vcf_list = ch_gatk_vcfs.collectFile( name: "vcf.list", newLine:true){ id,f -> "$f" } | collect

      gatk_mergevcfs(vcf_list,genome_dict)
    }
}

workflow {
// Prepare genome
  genome_fasta = Channel.fromPath(file(params.genome, checkIfExists:true)) | collect
  genome_index = bwa_index(genome_fasta) | collect
  genome_fai = faidx(genome_fasta) | collect
  genome_dict = gatk4_createsequencedict(genome_fasta) | collect


// Preprocess data
  ch_input_sample = extract_csv(file(params.samples, checkIfExists: true))

  ch_input_sample | qc

  ch_prep_reads = ch_input_sample | preprocess

  ch_mapped_marked_bams = gatk_map(ch_prep_reads,genome_fasta,genome_index, genome_dict)

  ch_mapped_marked_bais = ch_mapped_marked_bams | sidx

  ch_bbai = ch_mapped_marked_bams.join(ch_mapped_marked_bais)

  ch_bbai | bam_qc

// Parse caller param
  callerlist = params.callers?.split(',') as List

// Do the actual variant calling
  call_variants(callerlist,ch_bbai,genome_fasta,genome_fai,genome_dict)



}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        return [flowcell:fields[2],lane:fields[3]]
    } else if (fields.size() == 5) {
        return [flowcell:fields[0],lane:fields[1]]
    } else {
      // Not standard Illumina format
      // This print command is too spammy. TODO: Find a way to print just once at startup.
//      println "Sequence identifier does not conform to Illumina standard. Flowcell and lane will be set to dummy values"
      return [flowcell:"None",lane:1]
    }
}

def resolve_path(pathstring){
  if(pathstring =~ /^\//){
    pathstring
  } else {
    "${params.base_path}/${pathstring}"
  }
}

def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
    .map{ row -> 
      def meta = [:]
      meta.sample = row.sample

      def fastq_1     = file(resolve_path(row.fastq_1), checkIfExists: true)
      def fastq_2     = row.fastq_2 ? file(resolve_path(row.fastq_2), checkIfExists: true) : null
      meta.single_end = row.fastq_2 ? false : true

      def fclane = flowcellLaneFromFastq(fastq_1)
      meta.flowcell = fclane.flowcell
      meta.lane = fclane.lane  

      reads = [fastq_1,fastq_2]
      reads.removeAll([null])

      [meta,reads]
    }
}
