process multiqc_fastqc {

  publishDir "$params.outdir/multiqc", mode: 'copy'

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

  publishDir "$params.outdir/multiqc", mode: 'copy'

  input:
    path '*.json', stageAs: "?/*"

  output:
    path "fastp_multiqc_report.html", emit: report

  script:
  """
  multiqc -n fastp_multiqc_report.html .
  """
}



process multiqc_bams {

  publishDir "$params.outdir/multiqc", mode: 'copy'

  input:
    path '*.stat', stageAs: "?/*"
    path '*.flagstat', stageAs: "?/*"
    path '*.idxstat', stageAs: "?/*"

  output:
    path "bam_multiqc_report.html", emit: report

  script:
  """
  multiqc -n bam_multiqc_report.html .
  """
}
