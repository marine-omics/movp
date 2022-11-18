process fastp {

    publishDir "$params.outdir/fastp", mode: 'copy'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
        path('*.json')           , emit: json
        path('*.html')           , emit: html
        path('*.log')            , emit: log

    script:

        def read_group  = "${meta.sample}.${meta.flowcell}.${meta.lane}"
        def prefix = read_group
        def args = task.ext.args ?: ''

        if (meta.single_end) {
            """
            cat $reads \\
            | fastp \\
            --stdin \\
            --stdout \\
            --in1 ${reads[0]} \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $args \\
            2> ${prefix}.fastp.log \\
            | gzip -c > ${prefix}.fastp.fastq.gz
            """
            } else {
            """
            fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_1.fastp.fastq.gz \\
            --out2 ${prefix}_2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --thread $task.cpus \\
            $args \\
            2> ${prefix}.fastp.log
            """
        }
}

