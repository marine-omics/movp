# Marine Omics Variant Pipeline


# Quick Start

1. Install [nextflow](https://www.nextflow.io/)
2. Run a quick test to make sure everything is installed properly
```bash

```
3. Create the sample csv file


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
