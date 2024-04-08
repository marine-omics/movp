
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

Now supplementing this with a second region

```bash
# This is to find an appropriate contig
while read c;do cv=$(bcftools view 2785.vcf.gz -s "1142212,1142213,1142226,1142261,1142264,1142336,1142349,1151881"  -r $c | bcftools filter -i 'COUNT(GT="het")>0' | bcftools view -H | wc -l); echo $c,$cv  ;done < <(head -n 1000 contig_lengths.txt | cut -f 1)
```

```bash
for s in 1142212 1142213 1142226 1142261 1142264 1142336 1142349 1151881;do samtools view -O BAM ../${s}_mapped.bam QXJH01002620.1 > ${s}_mapped.bam; done

for f in *_mapped.bam;do samtools fastq $f | gzip > ${f%_mapped.bam}.fq.gz;done
```
