process multiqc_fastqc {

  publishDir "$params.outdir/multiqc"

  input:
    path '*.zip', stageAs: "?/*"

  output:
    path "fastqc_multiqc_report.html", emit: report

  script:
  """
  multiqc -n fastqc_multiqc_report.html .
  """
}


process multiqc_fastp {

  publishDir "$params.outdir/multiqc"

  input:
    path '*.json', stageAs: "?/*"

  output:
    path "fastp_multiqc_report.html", emit: report

  script:
  """
  multiqc -n fastp_multiqc_report.html .
  """
}
