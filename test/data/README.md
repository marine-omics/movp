
# Extract fastq from bam

```bash
for f in "DI-1-10.bam" "DI-1-1.bam";do
	sampleid=${f%.bam}
	samtools sort -n $f | samtools fastq -1 ${sampleid}_R1.fq -2 ${sampleid}_R2.fq -s /dev/null -
done

for fq in *.fq;do head -n 400 $fq | gzip > ${fq}.gz;done
```

And for SE data. Choose a few bams with a decent number of reads on the target locus

```bash
for f in *_mapped.bam;do cn=$(samtools view -c $f);echo "$f $cn";done | sort -n -k 2 | tail -n 8 > test_samples.txt
```

```bash
for f in $(cut -f 1 -d ' ' test_samples.txt);do
	samtools fastq $f | gzip > ${f%_mapped.bam}.fq.gz
done
```
