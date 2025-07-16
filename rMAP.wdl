version 1.0

workflow rMAP {
  input {
    Array[File] input_reads
    File adapters
    File reference_genome
    Boolean do_trimming = true
    Boolean do_quality_control = true
    Boolean do_assembly = true
    Boolean do_variant_calling = true
    Boolean do_annotation = true
    Boolean do_amr_profiling = true
    Boolean do_mlst = true
    Boolean do_pangenome = true
    Boolean do_phylogeny = true
    Boolean do_mge_analysis = true
    Boolean do_reporting = true
    Boolean do_blast = true
    Boolean use_local_blast = false
    File? local_blast_db
    String blast_db = "nt"
    Int blast_max_target_seqs = 250
    Float blast_evalue = 0.000001
    Int blast_min_contig_length = 300
    String? virulence_db = "vfdb"
    Int? virulence_min_cov = 60
    Float? virulence_min_id = 80.0
    String phylogeny_model = "-nt -gtr"
    String reference_type = "genbank"
    Int max_cpus = 16
    Int max_memory_gb = 32

    # Derived values explicitly computed in inputs section
    Int cpu_4 = if (max_cpus < 4) then max_cpus else 4
    Int cpu_8 = if (max_cpus < 8) then max_cpus else 8
    Int cpu_2 = if (max_cpus < 2) then max_cpus else 2
    Int mem_16 = if (max_memory_gb < 16) then max_memory_gb else 16
  }

  call CONFIGURATION {
    input:
      max_cpus = max_cpus,
      max_memory_gb = max_memory_gb
  }

  call TRIMMING {
    input:
      input_reads = input_reads,
      adapters = adapters,
      do_trimming = do_trimming,
      cpu = if 4 < max_cpus then 4 else max_cpus
  }

  call QUALITY_CONTROL {
    input:
      input_reads = if do_trimming then TRIMMING.trimmed_reads else input_reads,
      do_quality_control = do_quality_control,
      cpu = if 4 < max_cpus then 4 else max_cpus
  }

  call ASSEMBLY {
    input:
      input_reads = if do_trimming then TRIMMING.trimmed_reads else input_reads,
      assembler = "megahit",
      do_assembly = do_assembly,
      cpu = if 8 < max_cpus then 8 else max_cpus,
      memory_gb = if 16 < max_memory_gb then 16 else max_memory_gb
  }

  scatter (asm in ASSEMBLY.assembly_output) {
    call ASSEMBLY_STATS {
      input:
        assembly_file = asm,
        output_dir = ASSEMBLY.assembly_dir_out
    }
  }

  call ANNOTATION {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_annotation = do_annotation,
      cpu = cpu_4
  }

  scatter (sample in ANNOTATION.annotation_output) {
    call COPY_GFFS {
      input:
        gff_file = sample
    }
  }

  call MLST {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_mlst = do_mlst,
      cpu = cpu_2
  }

  call VARIANT_CALLING {
  input:
    input_reads = if do_trimming then TRIMMING.trimmed_reads else input_reads,
    reference_genome = reference_genome,
    do_variant_calling = do_variant_calling,
    reference_type = reference_type,
    cpu = cpu_4
  }

  call PANGENOME {
    input:
      annotation_input = flatten(COPY_GFFS.collected_gffs),
      do_pangenome = do_pangenome,
      cpu = cpu_4
  }

  call AMR_PROFILING {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_amr_profiling = do_amr_profiling,
      cpu = cpu_2
  }

  call MGE_ANALYSIS {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_mge_analysis = do_mge_analysis,
      cpu = cpu_4
  }

  call VIRULENCE_ANALYSIS {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      db_name = select_first([virulence_db]),
      min_coverage = select_first([virulence_min_cov]),
      min_identity = select_first([virulence_min_id]),
      cpu = cpu_2
  }

  call BLAST_ANALYSIS {
    input:
      contig_fastas = ASSEMBLY.assembly_output,
      blast_db = blast_db,
      max_target_seqs = blast_max_target_seqs,
      evalue = blast_evalue,
      min_contig_length = blast_min_contig_length,
      do_blast = do_blast
  }

  call PHYLOGENY {
    input:
      core_alignment = PANGENOME.core_alignment,
      accessory_alignment = PANGENOME.accessory_alignment,
      do_phylogeny = do_phylogeny,
      model = phylogeny_model,
      cpu = cpu_4
  }

  call REPORTING {
    input:
      mlst_output = select_first([MLST.combined_mlst]),
      amr_output = select_first([AMR_PROFILING.combined_amr]),
      phylogeny_output = PHYLOGENY.phylogeny_tree,
      variant_output = VARIANT_CALLING.vcf_files,
      blast_output = BLAST_ANALYSIS.blast_top5,
      plasmid_report = MGE_ANALYSIS.plasmid_report,
      virulence_report = VIRULENCE_ANALYSIS.combined_report,
      assembly_stats = ASSEMBLY_STATS.stats_txt
  }

  output {
    Array[File]? quality_reports = QUALITY_CONTROL.quality_reports
    Array[File]? trimmed_reads = TRIMMING.trimmed_reads
    Array[File]? assembly_output = ASSEMBLY.assembly_output
    Array[File]? assembly_stats = ASSEMBLY_STATS.stats_txt
    String? assembly_dir = ASSEMBLY.assembly_dir_out
    Array[File]? annotation_output = ANNOTATION.annotation_output
    Array[String]? annotation_dirs = ANNOTATION.annotation_dirs
    Array[File]? mlst_output = MLST.mlst_outputs
    File? combined_mlst = MLST.combined_mlst
    Array[File]? variant_output = VARIANT_CALLING.variant_output
    Array[String]? variant_dirs = VARIANT_CALLING.variant_dirs
    String? variants_dir = VARIANT_CALLING.variants_dir
    File? gene_presence_absence = PANGENOME.gene_presence_absence
    File? core_alignment = PANGENOME.core_alignment
    File? phylogeny_tree = PHYLOGENY.phylogeny_tree
    Array[File]? amr_output = AMR_PROFILING.amr_outputs
    File? combined_amr = AMR_PROFILING.combined_amr
    File? plasmid_report = MGE_ANALYSIS.plasmid_report
    Array[File]? plasmid_results = MGE_ANALYSIS.sample_reports
    Array[File]? blast_results = BLAST_ANALYSIS.blast_results
    Array[File]? blast_top5 = BLAST_ANALYSIS.blast_top5
    Array[File]? blast_logs = BLAST_ANALYSIS.blast_logs
    File? report_output = REPORTING.report_output
    Array[File]? virulence_reports = VIRULENCE_ANALYSIS.virulence_reports
    File? combined_virulence_report = VIRULENCE_ANALYSIS.combined_report
  }
}

task CONFIGURATION {
  input {
    Int max_cpus = 16
    Int max_memory_gb = 32
  }

  command {
    echo "System Configuration:"
    echo "Max CPUs: ~{max_cpus}"
    echo "Max Memory: ~{max_memory_gb} GB"
    echo "Workflow version: 1.1"
    date
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
  }

  output {
    String config_output = "Configuration complete"
  }
}

task TRIMMING {
  input {
    Array[File] input_reads
    File adapters
    Boolean do_trimming
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    if [ "~{do_trimming}" == "true" ]; then
      mkdir -p trimmed

      counter=0
      R1=""
      for file in ~{sep=' ' input_reads}; do
        if [ $((counter % 2)) -eq 0 ]; then
          R1="$file"
        else
          R2="$file"
          R1_base=$(basename "$R1")
          sample_name=$(echo "$R1_base" | sed -e 's/[._][Rr]1[._].*//' -e 's/[._]1[._].*//')

          echo "Processing sample: $sample_name (Files: $R1_base and $(basename "$R2"))"

          trimmomatic PE -threads ~{cpu} \
            "$R1" "$R2" \
            "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_1.unpair.fastq.gz" \
            "trimmed/${sample_name}_2.trim.fastq.gz" "trimmed/${sample_name}_2.unpair.fastq.gz" \
            ILLUMINACLIP:~{adapters}:2:30:10 \
            LEADING:20 \
            TRAILING:20 \
            SLIDINGWINDOW:4:20 \
            MINLEN:50

          for f in "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_2.trim.fastq.gz"; do
            if [ ! -f "$f" ]; then
              echo "Error: Failed to create $f" >&2
              exit 1
            fi
          done
        fi
        counter=$((counter + 1))
      done

      if [ $((counter % 2)) -ne 0 ]; then
        echo "ERROR: Odd number of input files - last file had no pair" >&2
        exit 1
      fi

      # Post-processing: Move and clean up
      echo "Moving trimmed files to main task directory..."
      mv trimmed ../trimmed

      echo "Cleaning up execution directory..."
      find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
    else
      echo "Trimming skipped by user request" > trimming_skipped.txt
    fi
  >>>

  runtime {
    docker: "staphb/trimmomatic:0.39"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: false
    preemptible: 2
  }

  output {
    Array[File] trimmed_reads = if do_trimming then glob("../trimmed/*_[12].trim.fastq.gz") else []
  }
}
task QUALITY_CONTROL {
  input {
    Array[File] input_reads
    Boolean do_quality_control
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    if [ "~{do_quality_control}" == "true" ]; then
      mkdir -p qc_reports
      echo "Running quality control on input reads..."
      fastqc -o qc_reports -t ~{cpu} ~{sep=' ' input_reads}

      # Run multiqc if available
      if which multiqc >/dev/null 2>&1; then
        echo "Running multiqc to aggregate reports..."
        multiqc qc_reports -o qc_reports
      fi

      # Verify reports were generated
      if [ $(ls qc_reports/*.{html,zip,txt} 2>/dev/null | wc -l) -eq 0 ]; then
        echo "ERROR: No QC reports generated!" >&2
        exit 1
      fi

      # Move reports to main task directory
      echo "Moving QC reports to main task directory..."
      mv qc_reports ../qc_reports

      # Clean up execution directory
      echo "Cleaning up execution directory..."
      find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
    else
      mkdir -p qc_reports
      echo "QC skipped by user request" > qc_reports/skipped.txt
      mv qc_reports ../qc_reports
    fi
  >>>

  runtime {
    docker: "staphb/fastqc:0.11.9"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
  }

  output {
    Array[File] quality_reports = if do_quality_control then glob("../qc_reports/*.{html,zip,txt}") else ["../qc_reports/skipped.txt"]
  }
}
task ASSEMBLY {
  input {
    Array[File] input_reads
    String assembler = "megahit"
    Boolean do_assembly = true
    String output_dir = "assembly"
    Int cpu = 4
    Int memory_gb = 16
  }

  command <<<
    set -euo pipefail

    if [ "~{do_assembly}" == "true" ]; then
      mkdir -p ~{output_dir}

      files=( ~{sep=' ' input_reads} )
      if [ $(( ${#files[@]} % 2 )) -ne 0 ]; then
          echo "ERROR: Odd number of input files for assembly" >&2
          exit 1
      fi

      for ((i=0; i<"${#files[@]}"; i+=2)); do
        R1="${files[i]}"
        R2="${files[i+1]}"
        sample_name=$(basename "$R1" | sed -E 's/_[R]?1[._].*//')
        outdir="~{output_dir}/megahit_${sample_name}"

        echo "Assembling $sample_name (R1: $R1, R2: $R2)"

        megahit \
          -1 "$R1" -2 "$R2" \
          -o "$outdir" \
          -t ~{cpu} \
          --memory 0.9 \
          --min-count 2 \
          --min-contig-len 500

        # Keep the original contigs.fa naming for consistency
        if [ -f "$outdir/final.contigs.fa" ]; then
          cp "$outdir/final.contigs.fa" "~{output_dir}/${sample_name}_contigs.fa"
        else
          echo "Warning: No contigs generated for $sample_name" >&2
          touch "~{output_dir}/${sample_name}_contigs.fa"
        fi
      done

      # Verify at least one contig file was created
      if [ $(find ~{output_dir} -name "*_contigs.fa" | wc -l) -eq 0 ]; then
        echo "ERROR: No contig files generated!" >&2
        exit 1
      fi

      # Move assembly directory to main task directory
      echo "Moving assembly results to main task directory..."
      mv ~{output_dir} ../~{output_dir}

      # Clean up execution directory
      echo "Cleaning up execution directory..."
      find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
    else
      mkdir -p ~{output_dir}
      echo "Assembly skipped by user request" > ~{output_dir}/skipped.txt
      mv ~{output_dir} ../~{output_dir}
    fi
  >>>

  runtime {
    docker: "quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 200 HDD"
    preemptible: 2
    continueOnReturnCode: true
  }

  output {
    Array[File] assembly_output = if do_assembly then glob("../~{output_dir}/*_contigs.fa") else []
    String assembly_dir_out = "../~{output_dir}"
  }
}

task ASSEMBLY_STATS {
  input {
    File assembly_file
    String output_dir
  }

  command <<<
    set -euo pipefail
    mkdir -p ~{output_dir}

    sample_name=$(basename "~{assembly_file}" .fa)
    echo "Calculating stats for: ~{assembly_file}"

    # Generate detailed stats (original full output)
    assembly-stats -t "~{assembly_file}" > "~{output_dir}/${sample_name}-stats.tab"

    echo "Assembly Statistics for ${sample_name}" > "~{output_dir}/${sample_name}-stats.txt"
    echo "=================================" >> "~{output_dir}/${sample_name}-stats.txt"
    assembly-stats "~{assembly_file}" >> "~{output_dir}/${sample_name}-stats.txt"
    echo "" >> "~{output_dir}/${sample_name}-stats.txt"
    echo "Contig Length Distribution (N50, L50, etc.):" >> "~{output_dir}/${sample_name}-stats.txt"
    echo "-------------------------------------------" >> "~{output_dir}/${sample_name}-stats.txt"
    assembly-stats --ng50 "~{assembly_file}" >> "~{output_dir}/${sample_name}-stats.txt"

    # Copy all stats files to parent directory (call-ASSEMBLY_STATS)
    find . -type f \( -name "*-stats.txt" -o -name "*-stats.tab" \) -exec cp -v {} .. \;

    # Safety check: Verify files exist in parent before deleting shards
    if ! ls ../*-stats.* >/dev/null 2>&1; then
      echo "ERROR: No stats files copied to parent directory!" >&2
      exit 1
    fi

    # Delete shard directories (only if copy succeeded)
    find . -type d -name "shard-*" -exec rm -rf {} +
  >>>

  runtime {
    docker: "quay.io/biocontainers/assembly-stats:1.0.1--h9948957_9"
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
  }

  output {
    File stats_tab = "~{output_dir}/${basename(assembly_file, '.fa')}-stats.tab"
    File stats_txt = "~{output_dir}/${basename(assembly_file, '.fa')}-stats.txt"
    Array[File] all_stats_files = glob("../*-stats.*")  # Captures all copied files
  }
}
task ANNOTATION {
  input {
    Array[File] assembly_output
    Boolean do_annotation = true
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    if [ "~{do_annotation}" == "true" ]; then
      mkdir -p annotation_results

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
        output_dir="annotation_results/${sample_name}"

        echo "Running PROKKA annotation for sample: $sample_name"

        prokka \
          --outdir "$output_dir" \
          --prefix "$sample_name" \
          --cpus ~{cpu} \
          --force \
          --compliant \
          --addgenes \
          "$asm_file" || {
            echo "PROKKA failed for $sample_name" >&2
            touch "$output_dir/$sample_name.gff"
          }
      done
    else
      echo "Annotation skipped by user request" > annotation_results/skipped.txt
    fi

    # Move results up one directory level and clean up
    mv annotation_results ../ && \
    cd ../ && \
    rm -rf execution || {
      echo "Warning: Failed to clean up execution directory" >&2
    }
  >>>

  runtime {
    docker: "staphb/prokka:1.14.6"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
    preemptible: 2
  }

  output {
    Array[File] annotation_output = if do_annotation then glob("annotation_results/*/*.gff") else []
    Array[String] annotation_dirs = if do_annotation then glob("annotation_results/*") else []
  }
}
task COPY_GFFS {
  input {
    File gff_file
  }

  command <<<
    set -euo pipefail
    mkdir -p ANNOTATED_GFFS
    cp "~{gff_file}" ANNOTATED_GFFS/
  >>>

  runtime {
    docker: "ubuntu:20.04"
    memory: "1 GB"
    cpu: 1
  }

  output {
    Array[File] collected_gffs = glob("ANNOTATED_GFFS/*")
  }
}

task MLST {
  input {
    Array[File] assembly_output
    Boolean do_mlst = true
    Int cpu = 2
  }

  command <<<
    set -euo pipefail

    if [ "~{do_mlst}" == "true" ]; then
      mkdir -p mlst_results

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
        output_file="mlst_results/${sample_name}_mlst.tsv"

        echo "Running MLST for sample: $sample_name"

        mlst "$asm_file" > "$output_file" || {
          echo "MLST failed for $sample_name" >&2
          echo -e "file\tscheme\tst\trep1\trep2\trep3\trep4\trep5\trep6\trep7" > "$output_file"
        }

        if [ -s "$output_file" ]; then
          awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
          mv "${output_file}.tmp" "$output_file"
        fi
      done

      if [ -n "$(ls -A mlst_results/*_mlst.tsv 2>/dev/null)" ]; then
        echo "Combining MLST results..."
        first_file=$(ls mlst_results/*_mlst.tsv | head -n1)
        head -n1 "$first_file" > mlst_results/combined_mlst.tsv
        for f in mlst_results/*_mlst.tsv; do
          tail -n +2 "$f" >> mlst_results/combined_mlst.tsv
        done
      else
        echo "No MLST results for any sample" > mlst_results/combined_mlst.tsv
      fi
    else
      echo "MLST skipped by user request" > mlst_results/skipped.txt
    fi

    # Move results up one directory level and clean up
    mv mlst_results ../ && \
    cd ../ && \
    rm -rf execution || {
      echo "Warning: Failed to clean up execution directory" >&2
    }
  >>>

  runtime {
    docker: "staphb/mlst:2.19.0"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
  }

  output {
    Array[File] mlst_outputs = if do_mlst then glob("mlst_results/*_mlst.tsv") else ["mlst_results/skipped.txt"]
    File? combined_mlst = if do_mlst then "mlst_results/combined_mlst.tsv" else "mlst_results/skipped.txt"
  }
}
task VARIANT_CALLING {
  input {
    Array[File] input_reads
    File reference_genome
    Boolean do_variant_calling = true
    String reference_type = "genbank"
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    if [ "~{do_variant_calling}" == "true" ]; then
      mkdir -p variants

      files=( ~{sep=' ' input_reads} )
      if [ $(( ${#files[@]} % 2 )) -ne 0 ]; then
          echo "ERROR: Odd number of input files" >&2
          exit 1
      fi

      for ((i=0; i<"${#files[@]}"; i+=2)); do
        R1="${files[i]}"
        R2="${files[i+1]}"
        SAMPLE_NAME=$(basename "$R1" | sed 's/_1.fastq.gz//; s/_1.trim.fastq.gz//')
        SAMPLE_DIR="variants/$SAMPLE_NAME"

        mkdir -p "$SAMPLE_DIR"
        cd "$SAMPLE_DIR"

        if [ "~{reference_type}" == "genbank" ]; then
          cp "~{reference_genome}" reference.gbk
          snippy \
            --cpus ~{cpu} \
            --ref reference.gbk \
            --R1 "$R1" \
            --R2 "$R2" \
            --outdir . \
            --prefix "$SAMPLE_NAME" \
            --force
        else
          cp "~{reference_genome}" reference.fasta
          samtools faidx reference.fasta
          snippy \
            --cpus ~{cpu} \
            --ref reference.fasta \
            --R1 "$R1" \
            --R2 "$R2" \
            --outdir . \
            --prefix "$SAMPLE_NAME" \
            --force
        fi

        cd ../..
      done

      # Verify variant files were created
      if [ $(find variants -name "*.snps.vcf" | wc -l) -eq 0 ]; then
        echo "ERROR: No VCF files generated!" >&2
        exit 1
      fi

      # Move variants directory to main task directory
      echo "Moving variant calling results to main task directory..."
      mv variants ../variants

      # Clean up execution directory
      echo "Cleaning up execution directory..."
      find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
    else
      mkdir -p variants
      echo "Variant calling skipped by user request" > variants/skipped.txt
      mv variants ../variants
    fi
  >>>

  runtime {
    docker: "quay.io/biocontainers/snippy:4.6.0--hdfd78af_1"
    memory: "16 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
    preemptible: 2
  }

  output {
    Array[File] vcf_files = if do_variant_calling then glob("../variants/*/*.snps.vcf") else []
    Array[File] variant_output = if do_variant_calling then glob("../variants/*/*") else []
    Array[String] variant_dirs = if do_variant_calling then glob("../variants/*") else []
    String variants_dir = "../variants"
  }
}
task PANGENOME {
  input {
    Array[File] annotation_input
    Boolean do_pangenome = true
    String output_prefix = "pangenome"
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

    if [ "~{do_pangenome}" == "true" ]; then
      mkdir -p gff_inputs
      for gff in ~{sep=' ' annotation_input}; do
        awk 'BEGIN {RS="\r\n"} {print}' "$gff" > "gff_inputs/$(basename "$gff")"
      done

      echo "=== EXECUTING ROARY ===" >&2
      roary -f ~{output_prefix}_results \
            -p ~{cpu} \
            -e -n -v \
            -i 90 -cd 99 \
            gff_inputs/*.gff || {
              echo "Roary failed with exit code $?" >&2
              exit 1
            }

      echo "=== LOCATING OUTPUT FILES ===" >&2
      output_dir="~{output_prefix}_results"

      [ -f "$output_dir/core_gene_alignment.aln" ] || {
        echo ">empty_sequence" > "$output_dir/core_gene_alignment.aln"
      }

      [ -f "$output_dir/accessory_gene_alignment.aln" ] || {
        echo ">empty_sequence" > "$output_dir/accessory_gene_alignment.aln"
      }

      [ -f "$output_dir/gene_presence_absence.csv" ] || {
        echo "Gene,Annotation" > "$output_dir/gene_presence_absence.csv"
      }

      mkdir -p final_output
      cp "$output_dir"/core_gene_alignment.aln final_output/
      cp "$output_dir"/accessory_gene_alignment.aln final_output/
      cp "$output_dir"/gene_presence_absence.csv final_output/
      cp "$output_dir"/accessory_binary_genes.fa final_output/ 2>/dev/null || touch final_output/accessory_binary_genes.fa
      cp "$output_dir"/summary_statistics.txt final_output/ 2>/dev/null || touch final_output/summary_statistics.txt

      # Move results to parent directory and clean up
      mv final_output ../ && \
      mv ~{output_prefix}_results ../ && \
      cd ../ && \
      rm -rf execution || {
        echo "Warning: Failed to clean up execution directory" >&2
      }
    else
      echo "Pangenome analysis skipped by user request" > skipped.txt
    fi
  >>>

  runtime {
    docker: "quay.io/biocontainers/roary:3.13.0--pl526h516909a_0@sha256:c188ba49c32c1a7204c64d5d3a14d4861c3b75741796b72dd4951e1d6af76bfd"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: false
    preemptible: 2
  }

  output {
      File gene_presence_absence = if do_pangenome then "final_output/gene_presence_absence.csv" else "skipped.txt"
      File summary_statistics = if do_pangenome then "final_output/summary_statistics.txt" else "skipped.txt"
      File accessory_binary = if do_pangenome then "final_output/accessory_binary_genes.fa" else "skipped.txt"
      File core_alignment = if do_pangenome then "final_output/core_gene_alignment.aln" else "skipped.txt"
      File accessory_alignment = if do_pangenome then "final_output/accessory_gene_alignment.aln" else "skipped.txt"
  }
}
task PHYLOGENY {
  input {
    File core_alignment
    File accessory_alignment
    Boolean do_phylogeny = true
    String tree_prefix = "phylogeny"
    String model = "-nt -gtr"
    Int cpu = 4
    Int bootstrap_replicates = 100  # New parameter for bootstrap replicates
  }

  command <<<
    set -euo pipefail

    if [ "~{do_phylogeny}" == "true" ]; then
      mkdir -p phylogeny_results

      function build_tree {
        local alignment=$1
        local output_prefix=$2
        local bootstrap=~{bootstrap_replicates}

        # Check if alignment exists and is not empty
        if [ ! -s "$alignment" ]; then
          echo "WARNING: Alignment file $alignment is empty or missing" >&2
          return 1
        fi

        # Count sequences in alignment
        seq_count=$(grep -c '^>' "$alignment" || echo 0)
        if [ "$seq_count" -lt 4 ]; then
          echo "WARNING: Insufficient sequences ($seq_count) in $alignment for phylogeny" >&2
          return 1
        fi

        echo "Building phylogeny for $output_prefix with $bootstrap bootstrap replicates..." >&2

        # Run FastTree with bootstrap support
        FastTree ~{model} \
          -gamma \
          -quiet \
          -boot $bootstrap \
          -log "phylogeny_results/${output_prefix}.log" \
          -out "phylogeny_results/${output_prefix}.nwk" \
          "$alignment" 2> "phylogeny_results/${output_prefix}.error.log" || {
            echo "FastTree failed for $output_prefix" >&2
            cat "phylogeny_results/${output_prefix}.error.log" >&2
            return 1
          }

        # Verify tree file was created
        if [ ! -s "phylogeny_results/${output_prefix}.nwk" ]; then
          echo "ERROR: Failed to generate tree file for $output_prefix" >&2
          return 1
        fi

        # Clean sample names in the tree file
        sed -i 's/\.contigs//g' "phylogeny_results/${output_prefix}.nwk"

        # Generate separate file with bootstrap values if needed
        if [ $bootstrap -gt 0 ]; then
          echo "Bootstrap support values included in the tree file" >&2
        fi
      }

      # Build core gene phylogeny
      echo "Building core gene phylogeny..." >&2
      build_tree "~{core_alignment}" "~{tree_prefix}_core" || {
        echo ">empty_sequence" > "phylogeny_results/~{tree_prefix}_core.nwk"
      }

      # Build accessory gene phylogeny
      echo "Building accessory gene phylogeny..." >&2
      build_tree "~{accessory_alignment}" "~{tree_prefix}_accessory" || {
        echo ">empty_sequence" > "phylogeny_results/~{tree_prefix}_accessory.nwk"
      }

      # Move all output files
      mv phylogeny_results/* ../ && \
      cd ../ && \
      rm -rf execution || {
        echo "Warning: Failed to clean up execution directory" >&2
      }
    else
      echo "Phylogeny skipped by user request" > skipped.txt
    fi
  >>>

  runtime {
    docker: "staphb/fasttree:2.1.11"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: [0, 1]
  }

  output {
    File phylogeny_tree = if do_phylogeny then "~{tree_prefix}_core.nwk" else "skipped.txt"
    File? phylogeny_log = if do_phylogeny then "~{tree_prefix}_core.log" else "skipped.txt"
    File? accessory_phylogeny = if do_phylogeny then "~{tree_prefix}_accessory.nwk" else "skipped.txt"
    File? core_error_log = if do_phylogeny then "~{tree_prefix}_core.error.log" else "skipped.txt"
    File? accessory_error_log = if do_phylogeny then "~{tree_prefix}_accessory.error.log" else "skipped.txt"
  }
}
task AMR_PROFILING {
  input {
    Array[File] assembly_output
    Boolean do_amr_profiling = true
    String db = "resfinder"
    Int minid = 90
    Int mincov = 80
    Int cpu = 2
  }

  command <<<
    set -euo pipefail

    if [ "~{do_amr_profiling}" == "true" ]; then
      mkdir -p amr_results

      if ! abricate --list | grep -q "~{db}"; then
        echo "ERROR: Database ~{db} not found" >&2
        exit 1
      fi

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
        output_file="amr_results/${sample_name}_amr.tsv"

        echo "Running AMR profiling for sample: $sample_name using ~{db} database"

        abricate \
          --db ~{db} \
          --minid ~{minid} \
          --mincov ~{mincov} \
          --threads ~{cpu} \
          "$asm_file" > "$output_file" || {
            echo "AMR profiling failed for $sample_name" >&2
            echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "$output_file"
          }

        if [ -s "$output_file" ]; then
          awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
          mv "${output_file}.tmp" "$output_file"
        fi
      done

      if [ -n "$(ls -A amr_results/*_amr.tsv 2>/dev/null)" ]; then
        echo "Combining AMR results..."
        first_file=$(ls amr_results/*_amr.tsv | head -n1)
        head -n1 "$first_file" > amr_results/combined_amr.tsv
        for f in amr_results/*_amr.tsv; do
          tail -n +2 "$f" >> amr_results/combined_amr.tsv
        done
      else
        echo "No AMR genes detected in any sample" > amr_results/combined_amr.tsv
      fi
    else
      echo "AMR profiling skipped by user request" > amr_results/skipped.txt
    fi

    # Move results up one directory level and clean up
    mv amr_results ../ && \
    cd ../ && \
    rm -rf execution || {
      echo "Warning: Failed to clean up execution directory" >&2
    }
  >>>

  runtime {
    docker: "staphb/abricate:1.0.0"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 20 HDD"
    continueOnReturnCode: true
  }

  output {
    Array[File] amr_outputs = if do_amr_profiling then glob("amr_results/*_amr.tsv") else ["amr_results/skipped.txt"]
    File? combined_amr = if do_amr_profiling then "amr_results/combined_amr.tsv" else "amr_results/skipped.txt"
  }
}
task MGE_ANALYSIS {
  input {
    Array[File] assembly_output
    Boolean do_mge_analysis = true
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    if [ "~{do_mge_analysis}" == "true" ]; then
      mkdir -p plasmid_results

      if ! abricate --list | grep -q plasmidfinder; then
        abricate --setupdb --db plasmidfinder || {
          echo "ERROR: Failed to setup plasmidfinder database" >&2
          exit 1
        }
      fi

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\..*//')
        echo "Running plasmid detection for: $sample_name"

        abricate \
          --db plasmidfinder \
          --mincov 80 \
          --minid 90 \
          --threads ~{cpu} \
          --nopath \
          "$asm_file" > "plasmid_results/${sample_name}.tsv" || {
            echo "Plasmid detection failed for $sample_name" >&2
            echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "plasmid_results/${sample_name}.tsv"
          }

        if [ -s "plasmid_results/${sample_name}.tsv" ]; then
          awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "plasmid_results/${sample_name}.tsv" > "plasmid_results/${sample_name}.tmp"
          mv "plasmid_results/${sample_name}.tmp" "plasmid_results/${sample_name}.tsv"
        fi
      done

      if [ -n "$(ls -A plasmid_results/*.tsv 2>/dev/null)" ]; then
        echo "Combining plasmid results..."
        first_file=$(ls plasmid_results/*.tsv | head -n1)
        head -n1 "$first_file" > plasmid_results/combined.tsv
        for f in plasmid_results/*.tsv; do
          tail -n +2 "$f" >> plasmid_results/combined.tsv
        done
      else
        echo "No plasmids detected in any sample" > plasmid_results/combined.tsv
      fi
    else
      echo "MGE analysis skipped by user request" > plasmid_results/skipped.txt
    fi

    # Move results up one directory level and clean up
    mv plasmid_results ../ && \
    cd ../ && \
    rm -rf execution || {
      echo "Warning: Failed to clean up execution directory" >&2
    }
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
  }

  output {
    File plasmid_report = if do_mge_analysis then "plasmid_results/combined.tsv" else "plasmid_results/skipped.txt"
    Array[File] sample_reports = if do_mge_analysis then glob("plasmid_results/*.tsv") else ["plasmid_results/skipped.txt"]
  }
}
task VIRULENCE_ANALYSIS {
  input {
    Array[File] assembly_output
    String db_name = "vfdb"
    Int min_coverage = 60
    Float min_identity = 80.0
    Int cpu = 2
  }

  command <<<
    set -euo pipefail

    if ! abricate --list | grep -q "~{db_name}"; then
      abricate --setupdb --db "~{db_name}" || {
        echo "ERROR: Failed to setup ~{db_name} database" >&2
        exit 1
      }
    fi

    mkdir -p virulence_results

    for asm_file in ~{sep=' ' assembly_output}; do
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//')
      output_file="virulence_results/${sample_name}_virulence.tsv"

      echo "Running virulence analysis for sample: $sample_name"

      abricate \
        --db ~{db_name} \
        --mincov ~{min_coverage} \
        --minid ~{min_identity} \
        --nopath \
        --quiet \
        --threads ~{cpu} \
        "$asm_file" > "$output_file" 2> "virulence_results/${sample_name}.log" || {
          echo "Abricate failed for $sample_name" >&2
          echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "$output_file"
        }

      if [ -s "$output_file" ]; then
        awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
        mv "${output_file}.tmp" "$output_file"
      fi
    done

    if [ -n "$(ls -A virulence_results/*_virulence.tsv 2>/dev/null)" ]; then
      echo "Combining virulence results..."
      first_file=$(ls virulence_results/*_virulence.tsv | head -n1)
      head -n1 "$first_file" > virulence_results/combined_virulence.tsv
      for f in virulence_results/*_virulence.tsv; do
        tail -n +2 "$f" >> virulence_results/combined_virulence.tsv
      done
    else
      echo "No virulence factors detected in any sample" > virulence_results/combined_virulence.tsv
    fi

    # Move results to main task directory
    echo "Moving virulence results to main task directory..."
    mkdir -p ../virulence_results
    mv virulence_results/*.tsv ../virulence_results/
    mv virulence_results/*.log ../virulence_results/

    # Clean up execution directory
    echo "Cleaning up execution directory..."
    find . -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    cpu: cpu
    memory: "4G"
    disks: "local-disk 20 HDD"
    preemptible: 2
    continueOnReturnCode: true
  }

  output {
    Array[File] virulence_reports = glob("../virulence_results/*_virulence.tsv")
    File combined_report = "../virulence_results/combined_virulence.tsv"
    Array[File] log_files = glob("../virulence_results/*.log")
  }
}

task BLAST_ANALYSIS {
  input {
    Array[File] contig_fastas
    String blast_db = "nt"
    Float evalue = 0.000001
    Int max_target_seqs = 250
    Int min_contig_length = 300
    Boolean do_blast = true
    Int cpu = 4
    Int memory_gb = 16
    Int max_retries_per_sample = 5
    Int retry_delay_seconds = 30
  }

  command <<<
    set -euo pipefail

    function extract_sample_id {
      basename "$1" | sed -E 's/([^._]+).*/\1/'
    }

    function filter_contigs {
      local input_file="$1"
      local min_length="$2"
      local output_file="$3"

      awk -v min_len="$min_length" 'BEGIN {RS=">";FS="\n"} {
          if (NR>1 && length($2) >= min_len) {
              print ">"$0
          }
      }' "$input_file" > "$output_file"

      [ -s "$output_file" ] || echo ">empty_sequence" > "$output_file"
    }

    function run_blast_with_retry {
      local query="$1"
      local output_file="$2"
      local log_file="$3"
      local attempt=1
      local max_attempts=~{max_retries_per_sample}
      local delay=~{retry_delay_seconds}

      while [ $attempt -le $max_attempts ]; do
        echo "Attempt $attempt of $max_attempts for $(basename $query)" >> "$log_file"

        if blastn \
          -query "$query" \
          -db "~{blast_db}" \
          -remote \
          -task blastn \
          -word_size 28 \
          -reward 1 -penalty -2 \
          -gapopen 2 -gapextend 1 \
          -outfmt "6 std qlen slen stitle" \
          -out "$output_file" \
          -evalue ~{evalue} \
          -max_target_seqs ~{max_target_seqs} \
          2>> "$log_file"; then

          if [ -s "$output_file" ] && grep -q -v '^#' "$output_file"; then
            echo "BLAST succeeded on attempt $attempt" >> "$log_file"
            return 0
          else
            echo "Empty/invalid BLAST output on attempt $attempt" >> "$log_file"
          fi
        else
          echo "BLAST failed with code $? on attempt $attempt" >> "$log_file"
        fi

        sleep $delay
        attempt=$((attempt + 1))
        delay=$((delay * 2))
      done

      echo "Max retries ($max_attempts) exceeded for $(basename $query)" >> "$log_file"
      return 1
    }

    if [ "~{do_blast}" = true ]; then
      mkdir -p blast_results
      > sample_ids.txt

      for contig_file in ~{sep=' ' contig_fastas}; do
        sample_id=$(extract_sample_id "$contig_file")
        filtered_contig="blast_results/${sample_id}_filtered.fa"
        blast_output="blast_results/${sample_id}_blast.tsv"
        blast_log="blast_results/${sample_id}_blast.log"

        echo "$sample_id" >> sample_ids.txt
        filter_contigs "$contig_file" ~{min_contig_length} "$filtered_contig"

        echo "Starting BLAST analysis for sample: $sample_id at $(date)" > "$blast_log"

        if run_blast_with_retry "$filtered_contig" "$blast_output" "$blast_log"; then
          echo "BLAST completed successfully for $sample_id" >> "$blast_log"
        else
          echo "Creating empty BLAST results after failures" >> "$blast_log"
          echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle" > "$blast_output"
        fi

        (
          head -n 1 "$blast_output" 2>/dev/null | grep -v '^#' || \
          echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle"
          grep -v '^#' "$blast_output" | \
          sort -t $'\t' -k1,1 -k12,12nr | \
          awk -F '\t' '
            BEGIN {OFS="\t"; prev=""}
            $1 != prev {count=1; prev=$1}
            count <= 5 {print; count++}
          '
        ) > "blast_results/${sample_id}_top5.tsv"
      done
    else
      mkdir -p blast_results
      echo "Remote BLAST skipped by user request" > blast_results/skipped.txt
      touch sample_ids.txt
    fi
  >>>

  runtime {
    docker: "ncbi/blast:2.14.0"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    maxRetries: 3
    preemptible: 2
    internet: true
    timeout: "12 hours"
  }

  output {
    Array[File] blast_results = if do_blast then glob("blast_results/*_blast.tsv") else []
    Array[File] blast_top5 = if do_blast then glob("blast_results/*_top5.tsv") else []
    Array[File] blast_logs = if do_blast then glob("blast_results/*_blast.log") else []
    Array[String] sample_ids = read_lines("sample_ids.txt")
  }
}
task REPORTING {
  input {
    File mlst_output
    File amr_output
    File phylogeny_output
    Array[File]? variant_output
    Array[File]? blast_output
    File plasmid_report
    File virulence_report
    Array[File] assembly_stats
  }

  command <<<
    set -euo pipefail

    function format_table {
      awk -F '\t' '
        BEGIN {
          OFS="\t"
          max_width=50
        }
        {
          for (i=1; i<=NF; i++) {
            if (length($i) > max_width) {
              $i=substr($i,1,max_width-3) "..."
            }
          }
          print
        }
      '
    }

    echo "Microbial Analysis Report" > report.txt
    echo "Generated: $(date)" >> report.txt
    echo "======================================" >> report.txt

    echo -e "\n=== Assembly Statistics ===" >> report.txt
    for stats_file in ~{sep=' ' assembly_stats}; do
      sample_name=$(basename "$stats_file" | sed 's/-stats.txt//')
      echo -e "\nSample: $sample_name" >> report.txt
      cat "$stats_file" >> report.txt
    done

    echo -e "\n=== MLST Results ===" >> report.txt
    cat ~{mlst_output} | format_table >> report.txt

    echo -e "\n=== AMR Profile ===" >> report.txt
    cat ~{amr_output} | format_table >> report.txt

    echo -e "\n=== Virulence Factors ===" >> report.txt
    cat ~{virulence_report} | format_table >> report.txt

    echo -e "\n=== Plasmid Detection ===" >> report.txt
    cat ~{plasmid_report} | format_table >> report.txt

    echo -e "\n=== Phylogenetic Tree ===" >> report.txt
    cat ~{phylogeny_output} >> report.txt

    if [ -n "$(ls -A ~{blast_output} 2>/dev/null)" ]; then
      echo -e "\n=== Top BLAST Hits ===" >> report.txt
      for blast_file in ~{sep=' ' blast_output}; do
        sample_name=$(basename "$blast_file" | sed 's/_top5.tsv//')
        echo -e "\nSample: $sample_name" >> report.txt
        cat "$blast_file" | format_table >> report.txt
      done
    fi

    if [ -n "$(ls -A ~{variant_output} 2>/dev/null)" ]; then
      echo -e "\n=== Variant Calling Results ===" >> report.txt
      for vcf_file in ~{sep=' ' variant_output}; do
        sample_name=$(basename "$vcf_file" | sed 's/.snps.vcf//')
        echo -e "\nSample: $sample_name" >> report.txt
        grep -v '^#' "$vcf_file" | head -n 20 | format_table >> report.txt
        echo "(showing first 20 variants)" >> report.txt
      done
    fi
  >>>

  runtime {
    docker: "ubuntu:20.04"
    memory: "4 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
    continueOnReturnCode: true
  }

  output {
    File report_output = "report.txt"
  }
}
