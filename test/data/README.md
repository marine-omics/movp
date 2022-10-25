
# Extract fastq from bam
```bash
f="DI-1-10.bam"
sampleid=${f%.bam}
samtools collate -u -O $f | samtools fastq -1 ${sampleid}_R1.fq -2 ${sampleid}_R2.fq -

for fq in *.fq;do head -n 400 $fq | gzip > ${fq}.gz;done
```