nextflow.enable.dsl=2

params.genome='../gatk/GCA_003918875.1_ASM391887v1_genomic.fna'
params.samples='samples.csv'
params.outdir='test'

process fastqc {

 //   container = 'marineomics/movp'
 // container 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0'
//  container 'quay.io/biocontainers/fastqc:0.11.9--0'
//  container 'file:///fast/sw/containers/fastqc-0.11.9.sif'

  input:
    path x

  output:
    path "*.zip"


  script:
  """
  fastqc $x
  """
}

process multiqc {

 // container = 'marineomics/movp'

  publishDir "$params.outdir/multiqc"
  input:
    path '*.zip', stageAs: "?/*"

  output:
    path "*multiqc_report.html", emit: report

  script:
  """
  multiqc  .
  """
}


workflow {
  ch_input_sample = extract_csv(file(params.samples, checkIfExists: true))

//  ch_input_sample | view { it.trim() }
  
  ch_input_sample | fastqc | collect | multiqc
  //| view { it.trim() }

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
    .map{ row -> row.fastq }
//    .view{ row -> row.fastq_1 }
//    .view { row }
}