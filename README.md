# Cardiogen WES Pipeline v2.0

![Nextflow](https://img.shields.io/badge/nextflow-23.10+-brightgreen)
![Container](https://img.shields.io/badge/container-Docker%20%7C%20Singularity-orange)
![License](https://img.shields.io/badge/license-MIT-blue)

A reproducible, scalable Nextflow pipeline for Whole Exome Sequencing (WES) data analysis with clinical-grade variant calling.


## Requirements

- **Java** 11 or higher
- **Nextflow** 23.10 or higher
- **Docker** or **Singularity**
- **Resources**: Minimum 8GB RAM, 4 CPU cores

## Reference Genome

The pipeline expects a pre-indexed **hg38 reference genome** in:
```
${projectDir}/reference/hg38/
```

**Required files:**
```
./reference/hg38/
├── hg38.analysisSet.fa      # Reference FASTA
├── hg38.analysisSet.fa.fai  # FASTA index
├── hg38.analysisSet.fa.amb  # BWA index files
├── hg38.analysisSet.fa.ann
├── hg38.analysisSet.fa.bwt
├── hg38.analysisSet.fa.pac
└── hg38.analysisSet.fa.sa

## Installation

### 1. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 2. Clone the repository
```bash
git clone https://github.com/nicolaedrabcinski/cardiogen_v2.git
cd cardiogen_v2
```

### 3. Build Docker image
```bash
docker build -t cardiogen-wes-tools:latest -f Dockerfile.tools .
```

## Quick Start
```bash
nextflow run main.nf \
  --input_dir /path/to/fastq \
  --outdir ./results \
  --genome GRCh38
```



## Pipeline Overview
```
FASTQ files → QC → Alignment → Variant Calling → Annotation → Report
```
