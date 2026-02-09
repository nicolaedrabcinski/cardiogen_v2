### Cardiogen WES Pipeline v2.0
https://img.shields.io/badge/nextflow-23.10+-brightgreen
https://img.shields.io/badge/container-Docker%2520%257C%2520Singularity-orange
https://img.shields.io/badge/license-MIT-blue

A reproducible, scalable Nextflow pipeline for Whole Exome Sequencing (WES) data analysis with clinical-grade variant calling.

#### Quick Start
git clone https://github.com/nicolaedrabcinski/cardiogen_v2.git
cd cardiogen_v2

#### Run the pipeline
```bash```
nextflow run main.nf --input_dir /path/to/fastq --outdir ./results

#### Requirements
- Java 11+ (required for Nextflow)
- Nextflow 23.10+
- Docker

#### Build the docker image
```bash```
docker build -t cardiogen-wes-tools:latest -f Dockerfile.tools .
