# rMAP-docker
This Docker image provides a ready-to-use environment for rMAP, a bioinformatics pipeline for analyzing microbial genomic data &amp; profiling AMR, Mobilome &amp; Virulome. It includes all required tools &amp; dependencies, enabling reproducible, scalable analysis of NGS data in research &amp; public health settings, particularly for low-resource environments.

# rMAP: Rapid Microbial Analysis Pipeline

![rMAP Logo](logo.jpg)

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

## Docker image🐳

The Docker image is hosted on [DockerHub](https://hub.docker.com/r/gmboowa/rmap):

```bash

docker pull gmboowa/rmap:1.0

docker run -it --rm gmboowa/rmap:1.0

```

---

## Usage

```bash

Let’s say you have this '~/rmap_project' directory on your computer:

```
Basic docker run command

```
docker run -it --rm \

  -v ~/rmap_project:/data \
  
  gmboowa/rmap:1.0 \
  
  rmap --config /data/config.json
  
  

gmboowa/rmap:1.0: Your Docker image.

rmap --config /data/config.json: Runs rMAP using the configuration file inside the container.

Sample 'config.json' file

Here’s a minimal example of a config.json needed for rMAP (adjust paths & parameters as needed):

json

{
  "input_files": [
    "samples/sample1_R1.fastq.gz",
    "samples/sample1_R2.fastq.gz"
  ],
  "output_dir": "results/",
  "threads": 8,
  "reference_genome": "refs/bacteria_ref.fasta",
  "amr_database": "dbs/resfinder.fasta",
  "taxonomic_classification": true,
  "quality_control": true,
  "assembly": true,
  "variant_calling": true
}

**Note**:

Place the samples/, refs/, & dbs/ folders inside your local ~/data directory.

results/ will be created inside that directory too.

-v ~/data:/data: Mounts your local folder into the container.

```

### Required options:

- `-i/--input`     Path to input raw reads (.fastq or .fastq.gz)
- `-o/--output`    Path to output directory
- `-r/--reference` Reference genome in `.gbk` format

### Optional:

- `-q/--trim`      Trimming quality threshold (default: 27)
- `-a/--assembly`  Assembly tool: `shovill` or `megahit`
- `-m/--amr`       Perform AMR profiling
- `-vc/--varcall`  Perform variant calling
- `-p/--phylogeny` Perform phylogenetic tree construction
- `-s/--pangenome` Perform pangenome analysis

---

## Authors

- [Gerald Mboowa](https://github.com/gmboowa)
- [Ivan Sserwadda](https://github.com/GunzIvan28)
- [Stephen Kanyerezi](https://github.com/Kanyerezi30)

---

## Resources

- GitHub: [https://github.com/GunzIvan28/rMAP](https://github.com/GunzIvan28/rMAP)
- Issues: [https://github.com/GunzIvan28/rMAP/issues](https://github.com/GunzIvan28/rMAP/issues)

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

