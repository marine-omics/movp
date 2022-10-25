process fastqc {

  input:
    tuple val(meta), path(reads)
//    path x

  output:
    path "*.zip"


  script:
  """
  fastqc $reads
  """
}