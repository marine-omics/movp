nextflow.enable.dsl=2

// Overridden in nextflow.config for test profiles
params.base_path="${launchDir}"

include { fastqc } from './modules/fastqc.nf'
include { multiqc } from './modules/multiqc.nf'


workflow {
  ch_input_sample = extract_csv(file(params.samples, checkIfExists: true))

  ch_input_sample | view

  ch_input_sample | fastqc | collect | multiqc
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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

      reads = [fastq_1,fastq_2]
      reads.removeAll([null])

      [meta,reads]
    }
}
