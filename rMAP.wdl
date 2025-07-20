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
    String virulence_db = "vfdb"
    Int virulence_min_cov = 60
    Float virulence_min_id = 80.0
    String phylogeny_model = "-nt -gtr"
    String reference_type = "genbank"
    Int max_cpus = 16
    Int max_memory_gb = 32
    Int min_assembly_quality = 50
    Int min_read_length = 50
    Int min_mapping_quality = 20

    # Derived values
    Int cpu_4 = if (max_cpus < 4) then max_cpus else 4
    Int cpu_8 = if (max_cpus < 8) then max_cpus else 8
    Int cpu_2 = if (max_cpus < 2) then max_cpus else 2
    Int mem_16 = if (max_memory_gb < 16) then max_memory_gb else 16
  }

  meta {
    workflow_timeout: "168 hours"
    workflow_heartbeat_interval: "10 minutes"
    workflow_heartbeat_ttl: "30 minutes"
    allowNestedInputs: true
    maxRetries: 3
    continueOnReturnCode: [0, 1]
    author: "Bioinformatics Team"
    email: "support@bioinfo.org"
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
      cpu = cpu_4,
      min_length = min_read_length
  }

  call QUALITY_CONTROL {
    input:
      input_reads = if do_trimming then TRIMMING.trimmed_reads else input_reads,
      do_quality_control = do_quality_control,
      cpu = cpu_4
  }

  call ASSEMBLY {
    input:
      input_reads = if do_trimming then TRIMMING.trimmed_reads else input_reads,
      assembler = "megahit",
      do_assembly = do_assembly,
      cpu = cpu_8,
      memory_gb = mem_16,
      min_quality = min_assembly_quality
  }

  call ANNOTATION {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_annotation = do_annotation,
      cpu = cpu_4
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
      cpu = cpu_4,
      min_quality = min_mapping_quality
  }

  call PANGENOME {
    input:
      annotation_input = ANNOTATION.annotation_output,
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
      db_name = virulence_db,
      min_coverage = virulence_min_cov,
      min_identity = virulence_min_id,
      cpu = cpu_2
  }

  call BLAST_ANALYSIS {
    input:
      contig_fastas = ASSEMBLY.assembly_output,
      blast_db = blast_db,
      max_target_seqs = blast_max_target_seqs,
      evalue = blast_evalue,
      min_contig_length = blast_min_contig_length,
      do_blast = do_blast,
      cpu = cpu_4,
      memory_gb = mem_16
  }

  call CORE_PHYLOGENY {
    input:
      alignment = PANGENOME.core_alignment,
      do_phylogeny = do_phylogeny,
      model = phylogeny_model,
      cpu = cpu_4,
      tree_prefix = "core_genes"
  }

  call ACCESSORY_PHYLOGENY {
    input:
      alignment = PANGENOME.accessory_binary,
      do_phylogeny = do_phylogeny,
      model = phylogeny_model,
      cpu = cpu_4,
      tree_prefix = "accessory_genes"
  }

  call REPORTING {
    input:
      mlst_output = MLST.combined_mlst,
      amr_output = AMR_PROFILING.combined_amr,
      core_phylogeny_output = select_first([CORE_PHYLOGENY.phylogeny_tree, "default_core_tree.nwk"]),
      accessory_phylogeny_output = select_first([ACCESSORY_PHYLOGENY.phylogeny_tree, "default_accessory_tree.nwk"]),
      variant_output = VARIANT_CALLING.vcf_files,
      blast_output = BLAST_ANALYSIS.blast_top5,
      plasmid_report = MGE_ANALYSIS.plasmid_report,
      virulence_report = VIRULENCE_ANALYSIS.combined_report
  }

  output {
    Array[File] quality_reports = QUALITY_CONTROL.quality_reports
    Array[File] trimmed_reads = TRIMMING.trimmed_reads
    Array[File] assembly_output = ASSEMBLY.assembly_output
    String assembly_dir = ASSEMBLY.assembly_dir_out
    Array[File] annotation_output = ANNOTATION.annotation_output
    Array[String] annotation_dirs = ANNOTATION.annotation_dirs
    Array[File] mlst_output = MLST.mlst_outputs
    File combined_mlst = MLST.combined_mlst
    Array[File] variant_output = VARIANT_CALLING.variant_output
    Array[String] variant_dirs = VARIANT_CALLING.variant_dirs
    String variants_dir = VARIANT_CALLING.variants_dir
    File gene_presence_absence = PANGENOME.gene_presence_absence
    File core_alignment = PANGENOME.core_alignment
    File accessory_alignment = PANGENOME.accessory_binary
    File core_phylogeny_output = select_first([CORE_PHYLOGENY.phylogeny_tree, "default_core_tree.nwk"])
    File accessory_phylogeny_output = select_first([ACCESSORY_PHYLOGENY.phylogeny_tree, "default_accessory_tree.nwk"])
    File? core_phylogeny_tree = CORE_PHYLOGENY.phylogeny_tree
    File? accessory_phylogeny_tree = ACCESSORY_PHYLOGENY.phylogeny_tree
    Array[File] amr_output = AMR_PROFILING.amr_outputs
    File combined_amr = AMR_PROFILING.combined_amr
    File plasmid_report = MGE_ANALYSIS.plasmid_report
    Array[File] plasmid_results = MGE_ANALYSIS.sample_reports
    Array[File] blast_results = BLAST_ANALYSIS.blast_results
    Array[File] blast_top5 = BLAST_ANALYSIS.blast_top5
    Array[File] blast_logs = BLAST_ANALYSIS.blast_logs
    File? report_output = REPORTING.report_output
    Array[File] virulence_reports = VIRULENCE_ANALYSIS.virulence_reports
    File combined_virulence_report = VIRULENCE_ANALYSIS.combined_report
  }
}


task CONFIGURATION {
  input {
    Int max_cpus = 16
    Int max_memory_gb = 32
  }

  command {
    echo "System Configuration:" > config.log
    echo "Max CPUs: ~{max_cpus}" >> config.log
    echo "Max Memory: ~{max_memory_gb} GB" >> config.log
    echo "Workflow version: 2.0" >> config.log
    echo "Start time: $(date)" >> config.log
    echo "Hostname: $(hostname)" >> config.log
    echo "Kernel: $(uname -a)" >> config.log
    cat config.log
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
    preemptible: 2
  }

  output {
    String config_output = read_string("config.log")
  }
}

task TRIMMING {
  input {
    Array[File]+ input_reads
    File adapters
    Boolean do_trimming
    Int cpu = 4
    Int min_length = 50
  }

  command <<<
    set -euo pipefail

    # Debugging: Log input files
    echo "Starting trimming process at $(date)" > trimming.log
    echo "Input files received:" >> trimming.log
    for f in ~{sep=' ' input_reads}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> trimming.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> trimming.log
    done

    if [ "~{do_trimming}" == "true" ]; then
      mkdir -p trimmed
      echo "Trimming parameters:" >> trimming.log
      echo "- CPU: ~{cpu}" >> trimming.log
      echo "- Min length: ~{min_length}" >> trimming.log

      counter=0
      R1=""
      for file in ~{sep=' ' input_reads}; do
        if [ $((counter % 2)) -eq 0 ]; then
          R1="$file"
        else
          R2="$file"
          R1_base=$(basename "$R1")
          sample_name=$(echo "$R1_base" | sed -e 's/[._][Rr]1[._].*//' -e 's/[._]1[._].*//')

          echo "Processing sample: $sample_name (Files: $R1_base and $(basename "$R2"))" >> trimming.log

          trimmomatic PE -threads ~{cpu} \
            "$R1" "$R2" \
            "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_1.unpair.fastq.gz" \
            "trimmed/${sample_name}_2.trim.fastq.gz" "trimmed/${sample_name}_2.unpair.fastq.gz" \
            ILLUMINACLIP:~{adapters}:2:30:10:8:true \
            LEADING:20 \
            TRAILING:20 \
            SLIDINGWINDOW:4:20 \
            MINLEN:~{min_length} 2>> trimming.log

          # Validate output files
          for f in "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_2.trim.fastq.gz"; do
            if [ ! -f "$f" ] || [ ! -s "$f" ]; then
              echo "Error: Failed to create valid output file $f" >> trimming.log
              exit 1
            fi
            echo "Created output file: $f ($(wc -c < "$f") bytes)" >> trimming.log
          done
        fi
        counter=$((counter + 1))
      done

      if [ $((counter % 2)) -ne 0 ]; then
        echo "ERROR: Odd number of input files - last file had no pair" >> trimming.log
        exit 1
      fi

      echo "Trimming completed successfully at $(date)" >> trimming.log
      echo "Output files created:" >> trimming.log
      ls -lh trimmed/* >> trimming.log
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
    timeout: "6 hours"
  }

  output {
    Array[File] trimmed_reads = if do_trimming then glob("trimmed/*_[12].trim.fastq.gz") else []
    File? trimming_log = if do_trimming then "trimming.log" else "trimming_skipped.txt"
    File? trimmed_files_list = if do_trimming then "trimmed/file_list.txt" else "trimming_skipped.txt"
  }
}

task QUALITY_CONTROL {
  input {
    Array[File]+ input_reads
    Boolean do_quality_control
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input files
    echo "Starting quality control at $(date)" > qc.log
    echo "Input files verification:" >> qc.log
    for f in ~{sep=' ' input_reads}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> qc.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> qc.log
    done

    if [ "~{do_quality_control}" == "true" ]; then
      mkdir -p qc_reports
      echo "Quality control parameters:" >> qc.log
      echo "- CPU: ~{cpu}" >> qc.log

      echo "Running FastQC..." >> qc.log
      fastqc -o qc_reports -t ~{cpu} ~{sep=' ' input_reads} 2>> qc.log || {
        echo "FastQC failed with exit code $?" >> qc.log
        exit 1
      }

      echo "Running MultiQC..." >> qc.log
      multiqc qc_reports -o qc_reports --force 2>> qc.log || {
        echo "MultiQC failed with exit code $?" >> qc.log
        exit 1
      }

      # Verify reports
      if [ $(ls qc_reports/*.{html,zip,txt} 2>/dev/null | wc -l) -eq 0 ]; then
        echo "ERROR: No QC reports generated!" >> qc.log
        exit 1
      fi

      echo "Quality control completed at $(date)" >> qc.log
      echo "Output files created:" >> qc.log
      ls -lh qc_reports/* >> qc.log
    else
      mkdir -p qc_reports
      echo "QC skipped by user request" > qc_reports/skipped.txt
    fi
  >>>

  runtime {
    docker: "staphb/fastqc:0.11.9"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
    timeout: "4 hours"
  }

  output {
    Array[File] quality_reports = if do_quality_control then glob("qc_reports/*.{html,zip,txt}") else glob("qc_reports/skipped.txt")
    File? qc_log = if do_quality_control then "qc.log" else "qc_reports/skipped.txt"
    File? qc_files_list = if do_quality_control then "qc_reports/file_list.txt" else "qc_reports/skipped.txt"
  }
}

task ASSEMBLY {
  input {
    Array[File]+ input_reads
    String assembler = "megahit"
    Boolean do_assembly = true
    String output_dir = "assembly"
    Int cpu = 4
    Int memory_gb = 16
    Int min_quality = 50
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input files
    echo "Starting assembly at $(date)" > assembly.log
    echo "Input files verification:" >> assembly.log
    for f in ~{sep=' ' input_reads}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> assembly.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> assembly.log
    done

    if [ "~{do_assembly}" == "true" ]; then
      mkdir -p ~{output_dir}
      echo "Assembly parameters:" >> assembly.log
      echo "- Assembler: ~{assembler}" >> assembly.log
      echo "- CPU: ~{cpu}" >> assembly.log
      echo "- Memory: ~{memory_gb} GB" >> assembly.log
      echo "- Min quality: ~{min_quality}" >> assembly.log

      files=( ~{sep=' ' input_reads} )
      if [ $(( ${#files[@]} % 2 )) -ne 0 ]; then
          echo "ERROR: Odd number of input files for assembly" >> assembly.log
          exit 1
      fi

      for ((i=0; i<"${#files[@]}"; i+=2)); do
        R1="${files[i]}"
        R2="${files[i+1]}"
        sample_name=$(basename "$R1" | sed -E 's/_[R]?1[._].*//')
        outdir="~{output_dir}/megahit_${sample_name}"

        echo "Assembling $sample_name (R1: $R1, R2: $R2)" >> assembly.log

        megahit \
          -1 "$R1" -2 "$R2" \
          -o "$outdir" \
          -t ~{cpu} \
          --memory 0.9 \
          --min-count 2 \
          --min-contig-len ~{min_quality} \
          --k-list 21,33,55,77,99,121 \
          --merge-level 20,0.98 \
          --prune-level 2 \
          --prune-depth 2 \
          2>> assembly.log

        if [ -f "$outdir/final.contigs.fa" ]; then
          cp "$outdir/final.contigs.fa" "~{output_dir}/${sample_name}_contigs.fa"
          contig_count=$(grep -c '^>' "~{output_dir}/${sample_name}_contigs.fa")
          echo "Generated contigs for $sample_name: $contig_count sequences" >> assembly.log
        else
          echo "Warning: No contigs generated for $sample_name" >> assembly.log
          touch "~{output_dir}/${sample_name}_contigs.fa"
        fi
      done

      if [ $(find ~{output_dir} -name "*_contigs.fa" | wc -l) -eq 0 ]; then
        echo "ERROR: No contig files generated!" >> assembly.log
        exit 1
      fi

      echo "Assembly completed at $(date)" >> assembly.log
      echo "Output files created:" >> assembly.log
      ls -lh ~{output_dir}/* >> assembly.log
    else
      mkdir -p ~{output_dir}
      echo "Assembly skipped by user request" > ~{output_dir}/skipped.txt
    fi
  >>>

  runtime {
    docker: "quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 200 HDD"
    preemptible: 2
    continueOnReturnCode: true
    timeout: "24 hours"
  }

  output {
    Array[File] assembly_output = if do_assembly then glob("~{output_dir}/*_contigs.fa") else []
    String assembly_dir_out = "~{output_dir}"
    File? assembly_log = if do_assembly then "assembly.log" else "~{output_dir}/skipped.txt"
    File? assembly_files_list = if do_assembly then "~{output_dir}/file_list.txt" else "~{output_dir}/skipped.txt"
  }
}

task ANNOTATION {
  input {
    Array[File]+ assembly_output
    Boolean do_annotation = true
    Int cpu = 4
  }

  command <<<
    set -euo pipefail
    mkdir -p annotation_output

    for asm_file in ~{sep=' ' assembly_output}; do
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
      mkdir -p "annotation_output/$sample_name"

      prokka \
        --outdir "annotation_output/$sample_name" \
        --prefix "$sample_name" \
        --cpus ~{cpu} \
        --force \
        "$asm_file" || {
          echo "Prokka failed, creating empty GFF" >&2
          echo "##gff-version 3" > "annotation_output/$sample_name/$sample_name.gff"
        }
    done
  >>>

  runtime {
    docker: "staphb/prokka:1.14.6"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    preemptible: 2
    timeout: "6 hours"
  }

  output {
    Array[File] annotation_output = if do_annotation then glob("annotation_output/*/*.gff") else []
    Array[String] annotation_dirs = if do_annotation then glob("annotation_output/*") else []
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

    # Create output directory
    mkdir -p "~{output_prefix}_results"

    echo "=== RUNNING ROARY ===" >&2
    echo "Input files: ~{sep=' ' annotation_input}" >&2

    # Run Roary with error handling
    roary -f "~{output_prefix}_results" \
          -p ~{cpu} \
          -e -n -v \
          -i 90 -cd 99 \
          ~{sep=' ' annotation_input} || {
            echo "WARNING: Roary completed with errors" >&2
            # Create minimal outputs if Roary fails
            echo ">empty_sequence" > "~{output_prefix}_results/core_gene_alignment.aln"
            echo ">empty_sequence" > "~{output_prefix}_results/accessory_binary_genes.fa"
            echo "Gene,Annotation" > "~{output_prefix}_results/gene_presence_absence.csv"
          }

    # Verify and standardize outputs
    echo "=== PROCESSING OUTPUTS ===" >&2
    mkdir -p final_output

    # Core alignment
    if [ -f "~{output_prefix}_results/core_gene_alignment.aln" ]; then
      cp "~{output_prefix}_results/core_gene_alignment.aln" final_output/
    else
      echo "WARNING: No core alignment found, creating empty file" >&2
      echo ">empty_sequence" > final_output/core_gene_alignment.aln
    fi

    # Accessory binary
    if [ -f "~{output_prefix}_results/accessory_binary_genes.fa" ]; then
      cp "~{output_prefix}_results/accessory_binary_genes.fa" final_output/
    else
      echo "WARNING: No accessory binary found, creating empty file" >&2
      echo ">empty_sequence" > final_output/accessory_binary_genes.fa
    fi

    # Other outputs
    cp "~{output_prefix}_results"/gene_presence_absence.csv final_output/ || \
      echo "Gene,Annotation" > final_output/gene_presence_absence.csv

    echo "=== FINAL OUTPUTS ===" >&2
    ls -l final_output >&2
  >>>

  runtime {
    docker: "quay.io/biocontainers/roary:3.13.0--pl526h516909a_0"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true  # Allow continuation even if Roary has warnings
  }

  output {
    File gene_presence_absence = "final_output/gene_presence_absence.csv"
    File accessory_binary = "final_output/accessory_binary_genes.fa"
    File core_alignment = "final_output/core_gene_alignment.aln"
    Boolean sufficient_samples = length(annotation_input) >= 3
  }
}
task CORE_PHYLOGENY {
  input {
    File alignment
    Boolean do_phylogeny = true
    String tree_prefix = "phylogeny"
    String model = "-nt -gtr"
    Int cpu = 4
    Int bootstrap_replicates = 100
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input file
    echo "Starting core genome phylogenetic analysis at $(date)" > phylogeny.log
    echo "Input file verification:" >> phylogeny.log
    if [ ! -f "~{alignment}" ]; then
      echo "ERROR: Alignment file not found" >> phylogeny.log
      exit 1
    fi
    echo "- ~{alignment} ($(wc -c < "~{alignment}") bytes)" >> phylogeny.log

    if [ "~{do_phylogeny}" == "true" ]; then
      echo "Phylogeny parameters:" >> phylogeny.log
      echo "- Model: ~{model}" >> phylogeny.log
      echo "- CPU: ~{cpu}" >> phylogeny.log
      echo "- Bootstrap replicates: ~{bootstrap_replicates}" >> phylogeny.log

      mkdir -p phylogeny_results

      # Validate input alignment
      if [ ! -s "~{alignment}" ]; then
        echo "ERROR: Core alignment file is empty or missing" >> phylogeny.log
        exit 1
      fi

      seq_count=$(grep -c '^>' "~{alignment}" || echo 0)
      if [ "$seq_count" -lt 4 ]; then
        echo "ERROR: Insufficient sequences ($seq_count) in alignment. Need at least 4." >> phylogeny.log
        exit 1
      fi

      # Run FastTree
      echo "Running FastTree for core genome..." >> phylogeny.log
      if ! FastTree ~{model} \
        -gamma \
        -quiet \
        -boot ~{bootstrap_replicates} \
        -log "phylogeny_results/core_~{tree_prefix}.log" \
        < "~{alignment}" > "phylogeny_results/core_~{tree_prefix}.nwk" 2>> phylogeny.log; then

        echo "FastTree failed. Error log:" >> phylogeny.log
        cat "phylogeny_results/core_error.log" >> phylogeny.log
        exit 1
      fi

      if [ ! -s "phylogeny_results/core_~{tree_prefix}.nwk" ]; then
        echo "ERROR: Core tree generation failed - empty output file" >> phylogeny.log
        exit 1
      fi

      # Create default tree if empty
      if [ ! -s "phylogeny_results/core_~{tree_prefix}.nwk" ]; then
        echo "Creating default core tree" >> phylogeny.log
        echo "(A,B);" > "phylogeny_results/core_~{tree_prefix}.nwk"
      fi

      echo "Core genome phylogenetic analysis completed at $(date)" >> phylogeny.log
      echo "Output files created:" >> phylogeny.log
      ls -lh phylogeny_results/* >> phylogeny.log
    else
      echo "Core phylogeny skipped by user request" > skipped.txt
      mkdir -p phylogeny_results
      echo "(A,B);" > "phylogeny_results/core_~{tree_prefix}.nwk"
    fi
  >>>

  runtime {
    docker: "staphb/fasttree:2.1.11"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: [0, 1]
    timeout: "12 hours"
  }

  output {
    File phylogeny_tree = if do_phylogeny then "phylogeny_results/core_~{tree_prefix}.nwk" else "phylogeny_results/core_~{tree_prefix}.nwk"
    File? phylogeny_log = if do_phylogeny then "phylogeny_results/core_~{tree_prefix}.log" else "skipped.txt"
    File? error_log = if do_phylogeny then "phylogeny_results/core_error.log" else "skipped.txt"
    File? execution_log = if do_phylogeny then "phylogeny.log" else "skipped.txt"
  }
}

task ACCESSORY_PHYLOGENY {
  input {
    File alignment
    Boolean do_phylogeny = true
    String tree_prefix = "phylogeny"
    String model = "-nt -gtr"
    Int cpu = 4
    Int bootstrap_replicates = 100
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input file
    echo "Starting accessory genome phylogenetic analysis at $(date)" > phylogeny.log
    echo "Input file verification:" >> phylogeny.log
    if [ ! -f "~{alignment}" ]; then
      echo "ERROR: Alignment file not found" >> phylogeny.log
      exit 1
    fi
    echo "- ~{alignment} ($(wc -c < "~{alignment}") bytes)" >> phylogeny.log

    if [ "~{do_phylogeny}" == "true" ]; then
      echo "Phylogeny parameters:" >> phylogeny.log
      echo "- Model: ~{model}" >> phylogeny.log
      echo "- CPU: ~{cpu}" >> phylogeny.log
      echo "- Bootstrap replicates: ~{bootstrap_replicates}" >> phylogeny.log

      mkdir -p phylogeny_results

      # Validate input alignment
      if [ ! -s "~{alignment}" ]; then
        echo "ERROR: Accessory alignment file is empty or missing" >> phylogeny.log
        exit 1
      fi

      seq_count=$(grep -c '^>' "~{alignment}" || echo 0)
      if [ "$seq_count" -lt 4 ]; then
        echo "ERROR: Insufficient sequences ($seq_count) in alignment. Need at least 4." >> phylogeny.log
        exit 1
      fi

      # Run FastTree
      echo "Running FastTree for accessory genome..." >> phylogeny.log
      if ! FastTree ~{model} \
        -gamma \
        -quiet \
        -boot ~{bootstrap_replicates} \
        -log "phylogeny_results/accessory_~{tree_prefix}.log" \
        < "~{alignment}" > "phylogeny_results/accessory_~{tree_prefix}.nwk" 2>> phylogeny.log; then

        echo "FastTree failed. Error log:" >> phylogeny.log
        cat "phylogeny_results/accessory_error.log" >> phylogeny.log
        exit 1
      fi

      if [ ! -s "phylogeny_results/accessory_~{tree_prefix}.nwk" ]; then
        echo "ERROR: Accessory tree generation failed - empty output file" >> phylogeny.log
        exit 1
      fi

      # Create default tree if empty
      if [ ! -s "phylogeny_results/accessory_~{tree_prefix}.nwk" ]; then
        echo "Creating default accessory tree" >> phylogeny.log
        echo "(A,B);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      fi

      echo "Accessory genome phylogenetic analysis completed at $(date)" >> phylogeny.log
      echo "Output files created:" >> phylogeny.log
      ls -lh phylogeny_results/* >> phylogeny.log
    else
      echo "Accessory phylogeny skipped by user request" > skipped.txt
      mkdir -p phylogeny_results
      echo "(A,B);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
    fi
  >>>

  runtime {
    docker: "staphb/fasttree:2.1.11"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: [0, 1]
    timeout: "12 hours"
  }

  output {
    File phylogeny_tree = if do_phylogeny then "phylogeny_results/accessory_~{tree_prefix}.nwk" else "phylogeny_results/accessory_~{tree_prefix}.nwk"
    File? phylogeny_log = if do_phylogeny then "phylogeny_results/accessory_~{tree_prefix}.log" else "skipped.txt"
    File? error_log = if do_phylogeny then "phylogeny_results/accessory_error.log" else "skipped.txt"
    File? execution_log = if do_phylogeny then "phylogeny.log" else "skipped.txt"
  }
}

task MLST {
  input {
    Array[File]+ assembly_output
    Boolean do_mlst = true
    Int cpu = 2
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input files
    echo "Starting MLST analysis at $(date)" > mlst.log
    echo "Input files verification:" >> mlst.log
    for f in ~{sep=' ' assembly_output}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> mlst.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> mlst.log
    done

    if [ "~{do_mlst}" == "true" ]; then
      echo "MLST parameters:" >> mlst.log
      echo "- CPU: ~{cpu}" >> mlst.log

      mkdir -p mlst_results

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
        output_file="mlst_results/${sample_name}_mlst.tsv"

        echo "Processing $sample_name" >> mlst.log
        mlst "$asm_file" > "$output_file" 2>> mlst.log || {
          echo "MLST failed for $sample_name" >> mlst.log
          echo -e "file\tscheme\tst\trep1\trep2\trep3\trep4\trep5\trep6\trep7" > "$output_file"
        }

        if [ -s "$output_file" ]; then
          awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
          mv "${output_file}.tmp" "$output_file"
        fi
      done

      if [ -n "$(ls -A mlst_results/*_mlst.tsv 2>/dev/null)" ]; then
        echo "Combining MLST results..." >> mlst.log
        first_file=$(ls mlst_results/*_mlst.tsv | head -n1)
        head -n1 "$first_file" > mlst_results/combined_mlst.tsv
        for f in mlst_results/*_mlst.tsv; do
          tail -n +2 "$f" >> mlst_results/combined_mlst.tsv
        done
      else
        echo "No MLST results for any sample" > mlst_results/combined_mlst.tsv
      fi

      echo "MLST completed at $(date)" >> mlst.log
      echo "Output files created:" >> mlst.log
      ls -lh mlst_results/* >> mlst.log
    else
      echo "MLST skipped by user request" > mlst_results/skipped.txt
    fi
  >>>

  runtime {
    docker: "staphb/mlst:2.19.0"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
    timeout: "6 hours"
  }

  output {
    Array[File] mlst_outputs = if do_mlst then glob("mlst_results/*_mlst.tsv") else ["mlst_results/skipped.txt"]
    File combined_mlst = if do_mlst then "mlst_results/combined_mlst.tsv" else "mlst_results/skipped.txt"
    File? mlst_log = if do_mlst then "mlst.log" else "mlst_results/skipped.txt"
  }
}

task VARIANT_CALLING {
  input {
    Array[File]+ input_reads
    File reference_genome
    Boolean do_variant_calling = true
    String reference_type = "genbank"
    Int cpu = 4
    Int min_quality = 20
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input files
    echo "Starting variant calling at $(date)" > variant.log
    echo "Input files verification:" >> variant.log
    for f in ~{sep=' ' input_reads}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> variant.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> variant.log
    done
    if [ ! -f "~{reference_genome}" ]; then
      echo "ERROR: Reference genome not found" >> variant.log
      exit 1
    fi
    echo "- Reference: ~{reference_genome} ($(wc -c < "~{reference_genome}") bytes)" >> variant.log

    if [ "~{do_variant_calling}" == "true" ]; then
      echo "Variant calling parameters:" >> variant.log
      echo "- Reference type: ~{reference_type}" >> variant.log
      echo "- CPU: ~{cpu}" >> variant.log
      echo "- Min quality: ~{min_quality}" >> variant.log

      mkdir -p variants

      files=( ~{sep=' ' input_reads} )
      if [ $(( ${#files[@]} % 2 )) -ne 0 ]; then
          echo "ERROR: Odd number of input files" >> variant.log
          exit 1
      fi

      for ((i=0; i<"${#files[@]}"; i+=2)); do
        R1="${files[i]}"
        R2="${files[i+1]}"
        SAMPLE_NAME=$(basename "$R1" | sed 's/_1.fastq.gz//; s/_1.trim.fastq.gz//')
        SAMPLE_DIR="variants/$SAMPLE_NAME"

        mkdir -p "$SAMPLE_DIR"
        cd "$SAMPLE_DIR"

        echo "Processing $SAMPLE_NAME" >> ../variant.log

        if [ "~{reference_type}" == "genbank" ]; then
          cp "~{reference_genome}" reference.gbk
          snippy \
            --cpus ~{cpu} \
            --minqual ~{min_quality} \
            --ref reference.gbk \
            --R1 "$R1" \
            --R2 "$R2" \
            --outdir . \
            --prefix "$SAMPLE_NAME" \
            --force 2>> ../variant.log
        else
          cp "~{reference_genome}" reference.fasta
          samtools faidx reference.fasta
          snippy \
            --cpus ~{cpu} \
            --minqual ~{min_quality} \
            --ref reference.fasta \
            --R1 "$R1" \
            --R2 "$R2" \
            --outdir . \
            --prefix "$SAMPLE_NAME" \
            --force 2>> ../variant.log
        fi

        cd ../..
      done

      if [ $(find variants -name "*.snps.vcf" | wc -l) -eq 0 ]; then
        echo "ERROR: No VCF files generated!" >> variant.log
        exit 1
      fi

      echo "Variant calling completed at $(date)" >> variant.log
      echo "Output files created:" >> variant.log
      find variants -type f -exec ls -lh {} \; >> variant.log
    else
      mkdir -p variants
      echo "Variant calling skipped by user request" > variants/skipped.txt
    fi
  >>>

  runtime {
    docker: "quay.io/biocontainers/snippy:4.6.0--hdfd78af_1"
    memory: "16 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
    preemptible: 2
    timeout: "12 hours"
  }

  output {
    Array[File] vcf_files = if do_variant_calling then glob("variants/*/*.snps.vcf") else []
    Array[File] variant_output = if do_variant_calling then glob("variants/*/*") else []
    Array[String] variant_dirs = if do_variant_calling then glob("variants/*") else []
    String variants_dir = "variants"
    File? variant_log = if do_variant_calling then "variant.log" else "variants/skipped.txt"
  }
}

task AMR_PROFILING {
  input {
    Array[File]+ assembly_output
    Boolean do_amr_profiling = true
    String db = "resfinder"
    Int minid = 90
    Int mincov = 80
    Int cpu = 2
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input files
    echo "Starting AMR profiling at $(date)" > amr.log
    echo "Input files verification:" >> amr.log
    for f in ~{sep=' ' assembly_output}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> amr.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> amr.log
    done

    if [ "~{do_amr_profiling}" == "true" ]; then
      echo "AMR profiling parameters:" >> amr.log
      echo "- Database: ~{db}" >> amr.log
      echo "- Min identity: ~{minid}" >> amr.log
      echo "- Min coverage: ~{mincov}" >> amr.log
      echo "- CPU: ~{cpu}" >> amr.log

      mkdir -p amr_results

      if ! abricate --list | grep -q "~{db}"; then
        echo "ERROR: Database ~{db} not found" >> amr.log
        exit 1
      fi

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
        output_file="amr_results/${sample_name}_amr.tsv"

        echo "Processing $sample_name" >> amr.log
        abricate \
          --db ~{db} \
          --minid ~{minid} \
          --mincov ~{mincov} \
          --threads ~{cpu} \
          "$asm_file" > "$output_file" 2>> amr.log || {
            echo "AMR profiling failed for $sample_name" >> amr.log
            echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "$output_file"
          }

        if [ -s "$output_file" ]; then
          awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
          mv "${output_file}.tmp" "$output_file"
        fi
      done

      if [ -n "$(ls -A amr_results/*_amr.tsv 2>/dev/null)" ]; then
        echo "Combining AMR results..." >> amr.log
        first_file=$(ls amr_results/*_amr.tsv | head -n1)
        head -n1 "$first_file" > amr_results/combined_amr.tsv
        for f in amr_results/*_amr.tsv; do
          tail -n +2 "$f" >> amr_results/combined_amr.tsv
        done
      else
        echo "No AMR genes detected in any sample" > amr_results/combined_amr.tsv
      fi

      echo "AMR profiling completed at $(date)" >> amr.log
      echo "Output files created:" >> amr.log
      ls -lh amr_results/* >> amr.log
    else
      echo "AMR profiling skipped by user request" > amr_results/skipped.txt
    fi
  >>>

  runtime {
    docker: "staphb/abricate:1.0.0"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 20 HDD"
    continueOnReturnCode: true
    timeout: "6 hours"
  }

  output {
    Array[File] amr_outputs = if do_amr_profiling then glob("amr_results/*_amr.tsv") else ["amr_results/skipped.txt"]
    File combined_amr = if do_amr_profiling then "amr_results/combined_amr.tsv" else "amr_results/skipped.txt"
    File? amr_log = if do_amr_profiling then "amr.log" else "amr_results/skipped.txt"
  }
}

task MGE_ANALYSIS {
  input {
    Array[File]+ assembly_output
    Boolean do_mge_analysis = true
    Int cpu = 4
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input files
    echo "Starting MGE analysis at $(date)" > plasmid.log
    echo "Input files verification:" >> plasmid.log
    for f in ~{sep=' ' assembly_output}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> plasmid.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> plasmid.log
    done

    if [ "~{do_mge_analysis}" == "true" ]; then
      echo "MGE analysis parameters:" >> plasmid.log
      echo "- CPU: ~{cpu}" >> plasmid.log

      mkdir -p plasmid_results

      if ! abricate --list | grep -q plasmidfinder; then
        abricate --setupdb --db plasmidfinder >> plasmid.log 2>&1 || {
          echo "ERROR: Failed to setup plasmidfinder database" >> plasmid.log
          exit 1
        }
      fi

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\..*//')
        echo "Processing $sample_name" >> plasmid.log

        abricate \
          --db plasmidfinder \
          --mincov 80 \
          --minid 90 \
          --threads ~{cpu} \
          --nopath \
          "$asm_file" > "plasmid_results/${sample_name}.tsv" 2>> plasmid.log || {
            echo "Plasmid detection failed for $sample_name" >> plasmid.log
            echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "plasmid_results/${sample_name}.tsv"
          }

        if [ -s "plasmid_results/${sample_name}.tsv" ]; then
          awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "plasmid_results/${sample_name}.tsv" > "plasmid_results/${sample_name}.tmp"
          mv "plasmid_results/${sample_name}.tmp" "plasmid_results/${sample_name}.tsv"
        fi
      done

      if [ -n "$(ls -A plasmid_results/*.tsv 2>/dev/null)" ]; then
        echo "Combining plasmid results..." >> plasmid.log
        first_file=$(ls plasmid_results/*.tsv | head -n1)
        head -n1 "$first_file" > plasmid_results/combined.tsv
        for f in plasmid_results/*.tsv; do
          tail -n +2 "$f" >> plasmid_results/combined.tsv
        done
      else
        echo "No plasmids detected in any sample" > plasmid_results/combined.tsv
      fi

      echo "MGE analysis completed at $(date)" >> plasmid.log
      echo "Output files created:" >> plasmid.log
      ls -lh plasmid_results/* >> plasmid.log
    else
      echo "MGE analysis skipped by user request" > plasmid_results/skipped.txt
    fi
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
    timeout: "6 hours"
  }

  output {
    File plasmid_report = if do_mge_analysis then "plasmid_results/combined.tsv" else "plasmid_results/skipped.txt"
    Array[File] sample_reports = if do_mge_analysis then glob("plasmid_results/*.tsv") else ["plasmid_results/skipped.txt"]
    File? plasmid_log = if do_mge_analysis then "plasmid.log" else "plasmid_results/skipped.txt"
  }
}

task VIRULENCE_ANALYSIS {
  input {
    Array[File]+ assembly_output
    String db_name = "vfdb"
    Int min_coverage = 60
    Float min_identity = 80.0
    Int cpu = 2
  }

  command <<<
    set -euo pipefail

    # Debugging: Verify input files
    echo "Starting virulence analysis at $(date)" > virulence.log
    echo "Input files verification:" >> virulence.log
    for f in ~{sep=' ' assembly_output}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> virulence.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> virulence.log
    done

    echo "Virulence analysis parameters:" >> virulence.log
    echo "- Database: ~{db_name}" >> virulence.log
    echo "- Min coverage: ~{min_coverage}" >> virulence.log
    echo "- Min identity: ~{min_identity}" >> virulence.log
    echo "- CPU: ~{cpu}" >> virulence.log

    if ! abricate --list | grep -q "~{db_name}"; then
      abricate --setupdb --db "~{db_name}" >> virulence.log 2>&1 || {
        echo "ERROR: Failed to setup ~{db_name} database" >> virulence.log
        exit 1
      }
    fi

    mkdir -p virulence_results

    for asm_file in ~{sep=' ' assembly_output}; do
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//')
      output_file="virulence_results/${sample_name}_virulence.tsv"

      echo "Processing $sample_name" >> virulence.log
      abricate \
        --db ~{db_name} \
        --mincov ~{min_coverage} \
        --minid ~{min_identity} \
        --nopath \
        --quiet \
        --threads ~{cpu} \
        "$asm_file" > "$output_file" 2>> virulence.log || {
          echo "Abricate failed for $sample_name" >> virulence.log
          echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "$output_file"
        }

      if [ -s "$output_file" ]; then
        awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
        mv "${output_file}.tmp" "$output_file"
      fi
    done

    if [ -n "$(ls -A virulence_results/*_virulence.tsv 2>/dev/null)" ]; then
      echo "Combining virulence results..." >> virulence.log
      first_file=$(ls virulence_results/*_virulence.tsv | head -n1)
      head -n1 "$first_file" > virulence_results/combined_virulence.tsv
      for f in virulence_results/*_virulence.tsv; do
        tail -n +2 "$f" >> virulence_results/combined_virulence.tsv
      done
    else
      echo "No virulence factors detected in any sample" > virulence_results/combined_virulence.tsv
    fi

    echo "Virulence analysis completed at $(date)" >> virulence.log
    echo "Output files created:" >> virulence.log
    ls -lh virulence_results/* >> virulence.log

    # Move results to main task directory
    mkdir -p ../virulence_results
    mv virulence_results/*.tsv ../virulence_results/
    mv virulence.log ../virulence_results/
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    cpu: cpu
    memory: "4G"
    disks: "local-disk 20 HDD"
    preemptible: 2
    continueOnReturnCode: true
    timeout: "6 hours"
  }

  output {
    Array[File] virulence_reports = glob("../virulence_results/*_virulence.tsv")
    File combined_report = "../virulence_results/combined_virulence.tsv"
    File? virulence_log = "../virulence_results/virulence.log"
  }
}

task BLAST_ANALYSIS {
  input {
    Array[File]+ contig_fastas
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

    # Debugging: Verify input files
    echo "Starting BLAST analysis at $(date)" > blast.log
    echo "Input files verification:" >> blast.log
    for f in ~{sep=' ' contig_fastas}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> blast.log
        exit 1
      fi
      echo "- $f ($(wc -c < "$f") bytes)" >> blast.log
    done

    if [ "~{do_blast}" = true ]; then
      echo "BLAST parameters:" >> blast.log
      echo "- Database: ~{blast_db}" >> blast.log
      echo "- E-value: ~{evalue}" >> blast.log
      echo "- Max targets: ~{max_target_seqs}" >> blast.log
      echo "- Min contig length: ~{min_contig_length}" >> blast.log
      echo "- CPU: ~{cpu}" >> blast.log
      echo "- Memory: ~{memory_gb} GB" >> blast.log

      mkdir -p blast_results
      > sample_ids.txt

      for contig_file in ~{sep=' ' contig_fastas}; do
        sample_id=$(extract_sample_id "$contig_file")
        filtered_contig="blast_results/${sample_id}_filtered.fa"
        blast_output="blast_results/${sample_id}_blast.tsv"
        blast_log="blast_results/${sample_id}_blast.log"

        echo "$sample_id" >> sample_ids.txt
        filter_contigs "$contig_file" ~{min_contig_length} "$filtered_contig"

        echo "Processing $sample_id" >> blast.log
        run_blast_with_retry "$filtered_contig" "$blast_output" "$blast_log" || {
          echo "Creating empty BLAST results after failures" >> "$blast_log"
          echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle" > "$blast_output"
        }

        # Create top 5 hits file
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

      echo "BLAST analysis completed at $(date)" >> blast.log
      echo "Output files created:" >> blast.log
      ls -lh blast_results/* >> blast.log
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
    File? execution_log = if do_blast then "blast.log" else "blast_results/skipped.txt"
    Array[String] sample_ids = read_lines("sample_ids.txt")
  }
}

task REPORTING {
  input {
    File mlst_output
    File amr_output
    File core_phylogeny_output
    File accessory_phylogeny_output
    Array[File]? variant_output
    Array[File]? blast_output
    File plasmid_report
    File virulence_report
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

    function count_features {
      awk -F '\t' 'NR>1 {print $1}' "$1" | sort -u | wc -l
    }

    function count_sequences {
      grep -c '^>' "$1" 2>/dev/null || echo 0
    }

    # Debugging: Verify input files
    echo "Generating comprehensive report at $(date)" > report.log
    echo "Input files verification:" >> report.log
    echo "- MLST: ~{mlst_output} ($(wc -c < "~{mlst_output}") bytes)" >> report.log
    echo "- AMR: ~{amr_output} ($(wc -c < "~{amr_output}") bytes)" >> report.log
    echo "- Core phylogeny: ~{core_phylogeny_output} ($(wc -c < "~{core_phylogeny_output}") bytes)" >> report.log
    echo "- Accessory phylogeny: ~{accessory_phylogeny_output} ($(wc -c < "~{accessory_phylogeny_output}") bytes)" >> report.log
    echo "- Plasmid report: ~{plasmid_report} ($(wc -c < "~{plasmid_report}") bytes)" >> report.log
    echo "- Virulence report: ~{virulence_report} ($(wc -c < "~{virulence_report}") bytes)" >> report.log

    echo "Microbial Genomic Analysis Report" > report.txt
    echo "==========================================" >> report.txt
    echo "Generated: $(date)" >> report.txt
    echo "Analysis Version: 2.1" >> report.txt
    echo "==========================================" >> report.txt

    echo -e "\n=== SUMMARY ===" >> report.txt
    echo -e "MLST profiles: $(count_features ~{mlst_output})" >> report.txt
    echo -e "AMR genes: $(count_features ~{amr_output})" >> report.txt
    echo -e "Virulence factors: $(count_features ~{virulence_report})" >> report.txt
    echo -e "Plasmids detected: $(count_features ~{plasmid_report})" >> report.txt

    echo -e "\n=== MLST Results ===" >> report.txt
    head -n 20 ~{mlst_output} | format_table >> report.txt
    [ $(wc -l < ~{mlst_output}) -gt 20 ] && echo "... (showing first 20 profiles)" >> report.txt

    echo -e "\n=== AMR Profile ===" >> report.txt
    head -n 20 ~{amr_output} | format_table >> report.txt
    [ $(wc -l < ~{amr_output}) -gt 20 ] && echo "... (showing first 20 genes)" >> report.txt

    echo -e "\n=== Virulence Factors ===" >> report.txt
    head -n 20 ~{virulence_report} | format_table >> report.txt
    [ $(wc -l < ~{virulence_report}) -gt 20 ] && echo "... (showing first 20 factors)" >> report.txt

    echo -e "\n=== Plasmid Detection ===" >> report.txt
    head -n 20 ~{plasmid_report} | format_table >> report.txt
    [ $(wc -l < ~{plasmid_report}) -gt 20 ] && echo "... (showing first 20 plasmids)" >> report.txt

    echo -e "\n=== Core Genome Phylogeny ===" >> report.txt
    echo -e "\nCore Genome Alignment Statistics:" >> report.txt
    echo -e "Number of sequences: $(count_sequences ~{core_phylogeny_output})" >> report.txt
    echo -e "\nNewick Tree:" >> report.txt
    cat ~{core_phylogeny_output} >> report.txt

    echo -e "\n=== Accessory Genome Phylogeny ===" >> report.txt
    echo -e "\nAccessory Genome Alignment Statistics:" >> report.txt
    echo -e "Number of sequences: $(count_sequences ~{accessory_phylogeny_output})" >> report.txt
    echo -e "\nNewick Tree:" >> report.txt
    cat ~{accessory_phylogeny_output} >> report.txt

    if [ -n "$(ls -A ~{blast_output} 2>/dev/null)" ]; then
      echo -e "\n=== Top BLAST Hits ===" >> report.txt
      for blast_file in ~{sep=' ' blast_output}; do
        sample_name=$(basename "$blast_file" | sed 's/_top5.tsv//')
        echo -e "\nSample: $sample_name" >> report.txt
        head -n 5 "$blast_file" | format_table >> report.txt
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

    echo -e "\n=== Analysis Summary ===" >> report.txt
    echo "Analysis completed successfully at $(date)" >> report.txt
    echo "Report generation completed at $(date)" >> report.log
  >>>

  runtime {
    docker: "ubuntu:20.04"
    memory: "4 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
    continueOnReturnCode: true
    timeout: "2 hours"
  }

  output {
    File report_output = "report.txt"
    File? report_log = "report.log"
  }
}
