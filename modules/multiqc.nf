process multiqc {

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
