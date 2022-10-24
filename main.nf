nextflow.enable.dsl=2


include { fastqc } from './modules/fastqc.nf'
include { multiqc } from './modules/multiqc.nf'


workflow test {
  ch_input_sample = extract_csv(file(params.samples, checkIfExists: true))
  ch_input_sample | fastqc 
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
def resolve_path(pathstring){
  if(pathstring =~ /^\//){
    pathstring
  } else {
    "${launchDir}/${pathstring}"
  }
}

// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
    .map{ row -> resolve_path(row.fastq) }
//    .view{ row -> row.fastq_1 }
//    .view { row }
}
