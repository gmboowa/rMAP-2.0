# rMAP-WDL-Cromwell-Docker
This tool provides a ready-to-use environment for rMAP, a bioinformatics pipeline for analyzing microbial genomic data &amp; profiling AMR, Mobilome &amp; Virulome. It includes all required tools &amp; dependencies, enabling reproducible, scalable analysis of NGS data in research &amp; public health settings.


**rMAP** is a fully automated pipeline for profiling the resistome & other genomic features of ESKAPEE (*Enterococcus faecium*, *Staphylococcus aureus*, *Klebsiella pneumoniae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*, *Enterobacter* species & *Escherichia coli*) pathogens using whole-genome sequencing (WGS) paired-end reads.

---

## Features

- Quality control of raw sequencing reads
- Trimming of adapters & low-quality bases
- *De novo* assembly of genomes
- Genome annotation (Prokka)
- Detection of antimicrobial resistance genes (AMR)
- Variant calling & SNP analysis
- Phylogenetic tree construction (IQ-TREE)
- Pangenome analysis (Roary)
- Mobile Genetic Element (MGE) profiling
- SRA data download support

---

# rMAP-WDL-Cromwell-Docker

**Version:** 1.0  
**Pipeline Type:** WDL-based, Docker-enabled  
**Workflow Engine:** Cromwell

## Overview
**rMAP-WDL-Cromwell-Docker** is a containerized, modular & scalable workflow for microbial genomics that integrates trimming, quality control, *de novo* assembly, annotation, variant calling, snpeff-based variant annotation, MLST typing, AMR profiling, mobile genetic element analysis, pangenome analysis, phylogeny & reporting.

This pipeline is written in **Workflow Description Language (WDL)**, utilizes **Docker containers** for tool standardization & is designed to run on the **Cromwell execution engine**.

## Features
- Adapter trimming with Trimmomatic
- FastQC-based quality control
- Megahit-based genome assembly
- Prokka for genome annotation
- Snippy for variant calling
- SnpEff for variant annotation
- MLST profiling
- Roary for pangenome construction
- FastTree for phylogenetic inference
- Abricate for AMR profiling
- mob-suite for MGE detection
- Remote BLAST support to NCBI
- Auto-generated final report

---

## Requirements
- [**Cromwell** (v84 or newer)](https://github.com/broadinstitute/cromwell/releases)
- [**Docker**](https://www.docker.com/) installed & running
- Input data: Paired-end FASTQ files
- Reference genome (FASTA)
- Adapter sequence file (FASTA or TXT)

---

![rMAP Logo](logo.jpg)

## How to download & run

### Step 1: Clone the repository
```bash
git clone https://github.com/gmboowa/rMAP-WDL-Cromwell-Docker.git
cd rMAP-WDL-Cromwell-Docker
```

### Step 2: Prepare inputs

Edit the input JSON file (e.g., `inputs.json`) with paths to your:
- Paired-end reads
- [**Reference genome**](https://www.ncbi.nlm.nih.gov/datasets/genome/)
- Illumina Adapter file
- Flags for toggling steps (true/false)

## SnpEff database references

```bash
# _Escherichia coli_ (K-12 substrain, RefSeq GCA_000974405)
docker run --rm -it staphb/snpeff:latest snpeff databases | grep -i Escherichia_coli_k_12_gca_000974405

# _Staphylococcus aureus subsp. aureus_ (NCTC 8325)
docker run --rm -it staphb/snpeff:latest snpeff databases | grep -i Staphylococcus_aureus_subsp_aureus_nctc_8325

# _Klebsiella pneumoniae subsp. pneumoniae_ (HS11286, RefSeq GCA_000240185)
docker run --rm -it staphb/snpeff:latest snpeff databases | grep -i Klebsiella_pneumoniae_subsp_pneumoniae_hs11286_gca_000240185

# _Acinetobacter baumannii_ (ATCC 19606 CIP 70.34)
docker run --rm -it staphb/snpeff:latest snpeff databases | grep -i Acinetobacter_baumannii_atcc_19606_cip_70_34

# _Pseudomonas aeruginosa_ (PAO1, RefSeq GCA_000006765)
docker run --rm -it staphb/snpeff:latest snpeff databases | grep -i Pseudomonas_aeruginosa_pao1_gca_000006765

# _Enterococcus faecium_ (S447)
docker run --rm -it staphb/snpeff:latest snpeff databases | grep -i Enterococcus_faecium_s447

# _Enterobacter cloacae_ (RefSeq GCA_001276405)
docker run --rm -it staphb/snpeff:latest snpeff databases | grep -i Enterobacter_cloacae_subsp_cloacae_gca_001276405

```
---
### Step 3: Run the workflow


## Run the command

```bash
java -jar cromwell.jar run rMAP.wdl --inputs inputs.json
```

To run on a backend like SLURM or Google Cloud, configure `cromwell.conf` accordingly.

---
### Note on pangenome & phylogenetic tree construction

- Pangenome analysis (Roary): Requires at least 3 annotated genome assemblies (in GFF3 format) for meaningful core/accessory genome separation.

- Phylogenetic tree construction (FastTree): Minimum of 4 samples is recommended to create a useful & interpretable tree. With fewer genomes, tree resolution & branching may be trivial or misleading.
  
- After running this tool on *Klebsiella pneumoniae* & *Escherichia coli*, proceed to analyze the assembled genomes using [`kleborate_wf.wdl`](https://github.com/gmboowa/kleborate_wf.wdl) to enable comprehensive genomic characterization.



---

## Output structure
- `trimmed/` – trimmed FASTQ files
- `qc_reports/` – FastQC reports
- `assembly/` – final contigs from Megahit
- `annotation_results/` – Prokka annotations
- `mlst_results/` – MLST profiles
- `variants/` – VCFs from Snippy
- `annotated_vcfs/` – SnpEff output
- `amr_results/` – AMR gene matches
- `mge_results/` – MGE prediction
- `pangenome_results/` – Roary files
- `phylogeny_results/` – Newick trees
- `remote_blast_results/` – BLAST XML files
- `report.txt` – consolidated plain-text report

---

## Sample input JSON
```json
{
  "rMAP.input_reads": [
    "~/test_data/A55870_1.fastq.gz",
    "~/test_data/A55870_2.fastq.gz",
    "~/test_data/A55888_1.fastq.gz",
    "~/test_data/A55888_2.fastq.gz"
  ],
  "rMAP.adapters": "~/adapters.fa",
  "rMAP.reference_genome": "~/GCF_000016305.1_ASM1630v1_genomic.fa",
  "rMAP.snpeff_organism": "Klebsiella_pneumoniae_subsp_pneumoniae_hs11286_gca_000240185",
  "rMAP.do_trimming": true,
  "rMAP.do_quality_control": true,
  "rMAP.do_assembly": true,
  "rMAP.do_annotation": true,
  "rMAP.do_mlst": true,
  "rMAP.do_variant_calling": true,
  "rMAP.do_pangenome": false,
  "rMAP.do_phylogeny": true,
  "rMAP.do_amr_profiling": true,
  "rMAP.do_mge_analysis": true,
  "rMAP.do_reporting": true,
  "rMAP.assembler": "megahit",
  "rMAP.do_blast": true,
  "rMAP.blast_db": "nt",
  "rMAP.blast_max_target_seqs": 250,
  "rMAP.blast_evalue": 0.000001
}
```

---

## Tools used (with Docker images)
| Step                | Tool          | Docker Image                          |
|---------------------|---------------|----------------------------------------|
| Trimming            | Trimmomatic   | `staphb/trimmomatic:0.39`             |
| QC                  | FastQC        | `staphb/fastqc:0.11.9`                |
| Assembly            | Megahit       | `quay.io/biocontainers/megahit:1.2.9` |
| Annotation          | Prokka        | `staphb/prokka:1.14.6`                |
| Variant Calling     | Snippy        | `staphb/snippy:4.6.0`                 |
| Variant Annotation  | SnpEff        | `staphb/snpeff:latest`                |
| MLST                | MLST          | `staphb/mlst:2.19.0`                  |
| Pangenome           | Roary         | `staphb/roary:3.13.0`                 |
| Phylogeny           | FastTree      | `staphb/fasttree:2.1.11`              |
| AMR Profiling       | Abricate      | `staphb/abricate:1.0.0`               |
| MGE Analysis        | mob-suite     | `continuumio/miniconda3:latest`       |
| Remote BLAST        | BLAST+        | `ncbi/blast:2.14.0`                   |


---

## Authors & contributors

- [Gerald Mboowa](https://github.com/gmboowa)
- [Ivan Sserwadda](https://github.com/GunzIvan28)
- [Stephen Kanyerezi](https://github.com/Kanyerezi30)


## Resources

- GitHub: [https://github.com/GunzIvan28/rMAP](https://github.com/GunzIvan28/rMAP)


---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

