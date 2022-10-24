process fastqc {

  input:
    path x

  output:
    path "*.zip"


  script:
  """
  fastqc $x
  """
}