# Marine Omics Variant Pipeline


# Quick Start

1. Install [nextflow](https://www.nextflow.io/)
2. Run a quick test to make sure everything is installed properly
```bash

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
