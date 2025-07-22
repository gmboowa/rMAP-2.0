# rMAP-2.0
This tool provides a ready-to-use environment for rMAP, a bioinformatics pipeline for analyzing microbial genomic data &amp; profiling AMR, Mobilome &amp; Virulome. It includes all required tools &amp; dependencies, enabling reproducible, scalable analysis of NGS data in research &amp; public health settings.


**rMAP** is a fully automated pipeline for profiling the resistome & other genomic features of ESKAPEE (*Enterococcus faecium*, *Staphylococcus aureus*, *Klebsiella pneumoniae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*, *Enterobacter* species & *Escherichia coli*) pathogens using whole-genome sequencing (WGS) paired-end reads.

---
## Overview
**Version:** 1.0  
**Pipeline Type:** WDL-based, Docker-enabled  
**Workflow Engine:** Cromwell

**rMAP-WDL-Cromwell-Docker** is a containerized, modular & scalable workflow for microbial genomics that integrates trimming, quality control, *de novo* assembly, annotation, variant calling, snpeff-based variant annotation, MLST typing, AMR profiling, mobile genetic element analysis, pangenome analysis, phylogeny & reporting.

This pipeline is written in **Workflow Description Language (WDL)**, utilizes **Docker containers** for tool standardization & is designed to run on the **Cromwell execution engine**.

## Features
- Adapter trimming with Trimmomatic  
- FastQC-based quality control  
- Megahit-based genome assembly  
- Prokka for genome annotation  
- Snippy for variant calling  
- MLST profiling  
- Roary for pangenome construction  
- FastTree for phylogenetic inference  
- Abricate for AMR & Virulence profiling & MGE detection  
- Remote BLAST support to NCBI  
- Visualize phylogenetic trees 


---

## Requirements
- [**Cromwell** (v84 or newer)](https://github.com/broadinstitute/cromwell/releases)
- [**Docker**](https://www.docker.com/) installed & running
- Input data: Paired-end FASTQ files
- Reference genome (FASTA)
- Adapter sequence file (FASTA or TXT)

---

![rMAP-2.0](rMAP-2.0.png)

## How to download & run

### Step 1: Clone the repository
```bash
git clone https://github.com/gmboowa/rMAP-2.0.git
cd rMAP-2.0
```

### Step 2: Prepare inputs

Edit the input JSON file (e.g., `inputs.json`) with paths to your:
- Paired-end reads
- [**Reference genome**](https://www.ncbi.nlm.nih.gov/datasets/genome/)
- Illumina Adapter file
- Flags for toggling steps (true/false)

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
    "~/A55738_1.fastq.gz",
    "~/A55738_2.fastq.gz",
    "~/A55870_1.fastq.gz",
    "~/A55870_2.fastq.gz",
    "~/A55888_1.fastq.gz",
    "~/A55888_2.fastq.gz",
    "~/A55944_1.fastq.gz",
    "~/A55944_2.fastq.gz",
    "~/A55727_1.fastq.gz",
    "~/A55727_2.fastq.gz"
 ],
  "rMAP.adapters": "/Volumes/MBOOWA/test_data/adapters.fa",
  "rMAP.reference_genome": "/Volumes/MBOOWA/test_data/GCA_000016305.1.gbk",
  "rMAP.do_trimming": true,
  "rMAP.do_quality_control": true,
  "rMAP.do_assembly": true,
  "rMAP.do_variant_calling": true,
  "rMAP.do_annotation": true,
  "rMAP.do_amr_profiling": true,
  "rMAP.do_mlst": true,
  "rMAP.do_pangenome": true,
  "rMAP.do_phylogeny": true,
  "rMAP.do_mge_analysis": true,
  "rMAP.do_reporting": true,
  "rMAP.do_blast": true,
  "rMAP.use_local_blast": true,
  "rMAP.local_blast_db": "~/refseq/bacteria/eskapee_combined.fasta",
  "rMAP.local_amr_db": "~/abricate/db/resfinder_db/resfinder.fa",
  "rMAP.local_mge_db": "~/abricate/db/plasmidfinder/plasmidfinder.fa",
  "rMAP.local_virulence_db": "~/abricate/db/vfdb/vfdb.fa",
  "rMAP.blast_db": "nt",
  "rMAP.blast_max_target_seqs": 250,
  "rMAP.blast_evalue": 0.000001,
  "rMAP.blast_min_contig_length": 300,
  "rMAP.virulence_db": "vfdb",
  "rMAP.virulence_min_cov": 60,
  "rMAP.virulence_min_id": 80.0,
  "rMAP.phylogeny_model": "-nt -gtr",
  "rMAP.reference_type": "genbank",
  "rMAP.max_cpus": 8,
  "rMAP.max_memory_gb": 16
}

```

---

## Tools used (with Docker images)
| Step                | Tool          | Docker image                          |
|---------------------|---------------|----------------------------------------|
| Trimming            | Trimmomatic   | `staphb/trimmomatic:0.39`             |
| QC                  | FastQC        | `staphb/fastqc:0.11.9`                |
| Assembly            | Megahit       | `quay.io/biocontainers/megahit:1.2.9` |
| Annotation          | Prokka        | `staphb/prokka:1.14.6`                |
| Variant Calling     | Snippy        | `staphb/snippy:4.6.0`                 |
| MLST                | MLST          | `staphb/mlst:2.19.0`                  |
| Pangenome           | Roary         | `staphb/roary:3.13.0`                 |
| Phylogeny           | FastTree      | `staphb/fasttree:2.1.11`              |
| AMR Profiling       | Abricate      | `staphb/abricate:1.0.0`               |
| MGE Analysis        | Abricate      | `staphb/abricate:latest`              |
| Virulence Analysis  | Abricate      | `staphb/abricate:latest`              |
| Remote BLAST        | BLAST+        | `ncbi/blast:2.14.0`                   |



---


## Output directory hierarchy for rMAP 2.0 Pipeline

After successful execution of the `rMAP` pipeline using WDL + Cromwell + Docker, your output directory will contain subdirectories corresponding to each major analysis module. Below is the typical hierarchy:

```bash
rMAP_outputs/
├── call-CONFIGURATION/
├── call-TRIMMING/
├── call-QUALITY_CONTROL/
├── call-ASSEMBLY/
├── call-VARIANT_CALLING/
├── call-AMR_PROFILING/
├── call-MLST/
├── call-MGE_ANALYSIS/
├── call-VIRULENCE_ANALYSIS/
├── call-ANNOTATION/
├── call-BLAST_ANALYSIS/
├── call-PANGENOME/
├── call-ACCESSORY_PHYLOGENY/
├── call-CORE_PHYLOGENY/
├── call-TREE_VISUALIZATION/
```

Each `call-*` directory contains:
- `execution/` – Shell scripts & logs for the task.
- `stdout` / `stderr` – Standard output & error logs.
- `rc` – Return code for the task.
- Output files generated by the task (e.g., `.fasta`, `.vcf`, `.tsv`, `.json`, `.html`, etc.).

## Example Outputs by Module

| Module                   | Key output files                                                 |
|--------------------------|------------------------------------------------------------------|
| `TRIMMING`               | Trimmed FASTQ files  (`*.fastq.gz`)                              |
| `QUALITY_CONTROL`        | MultiQC reports, FastQC files (`*.zip`, `*.html`)                |
| `ASSEMBLY`               | Assembled contigs (`*.fasta`)                                    |
| `VARIANT_CALLING`        | VCF files (`*.vcf`)                                              |
| `AMR_PROFILING`          | Resistance profiles (`*.txt`, `*.tsv`)                           |
| `MLST`                   | MLST profiles  (`*.txt`, `*.tsv`)                                |
| `MGE_ANALYSIS`           | Mobile genetic element annotations  (`*.txt`, `*.tsv`)           |
| `VIRULENCE_ANALYSIS`     | Virulence gene predictions  (`*.txt`, `*.tsv`)                   |
| `ANNOTATION`             | Genomic feature annotations (`*.gff`, `*.gbk`)                   |
| `BLAST_ANALYSIS`         | Top BLAST hits (`*.tsv`, `*.xml`)                                |
| `PANGENOME`              | Roary outputs: `gene_presence_absence.csv`, `core_gene_alignment.aln` |
| `ACCESSORY_PHYLOGENY`    | Phylogenetic tree for accessory genes (`*.nwk`, `*.pdf`)         |
| `CORE_PHYLOGENY`         | Core genome tree and alignment files (`*.nwk`, `*.pdf`)          |

## Note on BLAST usage
If you are analyzing many samples, we recommend setting up a local BLAST nucleotide database specifically for ESKAPEE pathogens. This setup requires approximately 70 GB of disk space.
Please note that NCBI imposes usage limits on BLAST queries from a single IP address, which may affect performance or availability during high-throughput runs. A local database ensures speed, reproducibility & compliance with query limits.

## Note on MLST Schemas
If you are performing MLST typing across many samples, we recommend downloading & setting up the publicly available PubMLST schemes locally. This setup requires approximately 2 GB of disk space. A local installation ensures faster typing, avoids dependency on internet connectivity & supports reproducible & scalable analysis across multiple species.

## Indexing custom BLAST databases

Before running rMAP v2.0, if you intend to use local databases, make sure to index the custom nucleotide FASTA files used for resistance, plasmid & virulence factor detection using `makeblastdb`. This step is necessary to enable local BLAST searches during the workflow. For each database file, run the following command to generate the required BLAST index files:


```bash
makeblastdb -in resfinder.fa -dbtype nucl
makeblastdb -in plasmidfinder.fa -dbtype nucl
makeblastdb -in vfdb.fa -dbtype nucl
```

## Authors & contributors

- [Gerald Mboowa](https://github.com/gmboowa)
- [Ivan Sserwadda](https://github.com/GunzIvan28)
- [Stephen Kanyerezi](https://github.com/Kanyerezi30)


## Resources

- GitHub: [https://github.com/GunzIvan28/rMAP](https://github.com/GunzIvan28/rMAP)
- Publication: rMAP: the Rapid Microbial Analysis Pipeline for ESKAPE bacterial group whole-genome sequence data  
*Published in [Microbial Genomics](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000583)*
- The philosophy of **rMAP-WDL-Cromwell-Docker** is built on the foundation of already pre-existing tools. As a token of gratitude to the authors of those numerous tools.
---

## License

This project is licensed under the MIT License. 

## To report bugs, ask questions or seek help

The software developing team works round the clock to ensure the bugs within the tool are captured and fixed. For support or any inquiry: You can submit your query using the [Issue Tracker](https://github.com/gmboowa/rMAP-WDL-Cromwell-Docker/issues)

