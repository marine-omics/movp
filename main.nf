nextflow.enable.dsl=2

include { fastqc } from './modules/fastqc.nf'
include { multiqc_fastqc; multiqc_fastp; multiqc_bams } from './modules/multiqc.nf'
include { fastp } from './modules/fastp.nf'
include { fastq2ubam; markadapters; bwa_mem_gatk; gatk4_createsequencedict; gatk_mark_duplicates } from './modules/gatk.nf'
include { bwa_index } from './modules/bwa.nf'
include { sidx; faidx; flagstat; stat; idxstat } from './modules/samtools.nf'
include { freebayes } from './modules/freebayes.nf'


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

    ch_merge_bams.view()

    mapped_bams = bwa_mem_gatk(ch_merge_bams,genome_fasta,genome_index, genome_dict)
    mapped_marked_bams = gatk_mark_duplicates(mapped_bams)

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

  ch_bbai_collection = ch_mapped_marked_bams.join(ch_mapped_marked_bais)

  ch_bbai_collection | bam_qc

// Freebayes
  ch_bamcollection = ch_mapped_marked_bams.map{m,b -> b} | collect
  ch_baicollection = ch_mapped_marked_bais.map{m,b -> b} | collect
  freebayes(ch_bamcollection,ch_baicollection,genome_fasta,genome_fai,file(params.populations))

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
