# Marine Omics Variant Pipeline


# Quick Start

1. Install [nextflow](https://www.nextflow.io/)
2. Run a test to make sure everything is installed properly. The command below should work on a linux machine with singularity installed. If you are working from a mac or windows machine you will need to use docker. 
```bash
nextflow run marine-omics/movp -profile singularity,test -r main
```
3. Create the sample csv file

Example (single end reads)
```
sample,fastq
1,sample1.fastq.gz
2,sample2.fastq.gz
```

Example (paired-end reads)
```
sample,fastq
1,sample1_r1.fastq.gz,sample1_r2.fastq.gz
2,sample2_r1.fastq.gz,sample2_r2.fastq.gz
```

4. Choose a profile for your execution environment


```bash
nextflow run marine-omics/movp -profile docker,test -r main
```

Paths should either be given as absolute paths or relative to the launch directory (where you invoked the nextflow command)

# Docker and Singularity

All of the dependencies for `movp` are provided in a single docker container. 

To build the container on a system with docker installed

```bash
	docker build -t marineomics/movp .
```

On a systems with singularity you will need to use this image as a `.sif` file.  To convert

```bash
docker save marineomics/movp -o movp.tar
```

Copy to a system with singularity installed and then run

```bash
singularity build movp.sif docker-archive://movp.tar
```
