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
    Boolean use_local_blast = true
    File? local_blast_db
    File? local_amr_db
    File? local_mge_db
    File? local_virulence_db
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
    Int tree_image_width = 1200
    String tree_image_format = "png"
    Int tree_font_size = 8

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
      cpu = cpu_4,
      output_prefix = "rMAP_pangenome"
  }

  call AMR_PROFILING {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_amr_profiling = do_amr_profiling,
      local_db = local_amr_db,
      use_local_db = use_local_blast,
      cpu = cpu_2
  }

  call MGE_ANALYSIS {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_mge_analysis = do_mge_analysis,
      local_db = local_mge_db,
      use_local_db = use_local_blast,
      cpu = cpu_4
  }

  call VIRULENCE_ANALYSIS {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      db_name = virulence_db,
      local_db = local_virulence_db,
      use_local_db = use_local_blast,
      min_coverage = virulence_min_cov,
      min_identity = virulence_min_id,
      cpu = cpu_2
  }

  call BLAST_ANALYSIS {
    input:
      contig_fastas = ASSEMBLY.assembly_output,
      blast_db = blast_db,
      local_blast_db = local_blast_db,
      use_local_blast = use_local_blast,
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
      tree_prefix = "core_genes",
      bootstrap_replicates = 100
  }

  call ACCESSORY_PHYLOGENY {
    input:
      alignment = PANGENOME.accessory_binary,
      do_phylogeny = do_phylogeny,
      model = phylogeny_model,
      cpu = cpu_4,
      tree_prefix = "accessory_genes",
      bootstrap_replicates = 100
  }

  # Tree visualization section - fixed implementation
  Array[File] phylogeny_trees = [
    CORE_PHYLOGENY.phylogeny_tree,
    ACCESSORY_PHYLOGENY.phylogeny_tree
  ]

  scatter (tree in phylogeny_trees) {
    call TREE_VISUALIZATION as TREE_VISUALIZATION {
      input:
        input_tree = tree,
        width = tree_image_width,
        image_format = tree_image_format,
        font_size = tree_font_size
    }
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
    Array[File] amr_output = AMR_PROFILING.amr_outputs
    File combined_amr = AMR_PROFILING.combined_amr
    File plasmid_report = MGE_ANALYSIS.plasmid_report
    Array[File] plasmid_results = MGE_ANALYSIS.sample_reports
    Array[File] blast_results = BLAST_ANALYSIS.blast_results
    Array[File] blast_top5 = BLAST_ANALYSIS.blast_top5
    Array[File] blast_logs = BLAST_ANALYSIS.blast_logs
    Array[File] virulence_reports = VIRULENCE_ANALYSIS.virulence_reports
    File combined_virulence_report = VIRULENCE_ANALYSIS.combined_report
    Array[File?] tree_images = TREE_VISUALIZATION.final_image
    Array[File] tree_render_logs = TREE_VISUALIZATION.render_log
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

          # Run trimmomatic in a subshell to prevent failure from stopping the whole task
          (
            trimmomatic PE -threads ~{cpu} \
              "$R1" "$R2" \
              "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_1.unpair.fastq.gz" \
              "trimmed/${sample_name}_2.trim.fastq.gz" "trimmed/${sample_name}_2.unpair.fastq.gz" \
              ILLUMINACLIP:~{adapters}:2:30:10:8:true \
              LEADING:20 \
              TRAILING:20 \
              SLIDINGWINDOW:4:20 \
              MINLEN:~{min_length} 2>> trimming.log || {
                echo "WARNING: Trimmomatic failed for sample $sample_name" >> trimming.log
                # Create empty output files to allow pipeline to continue
                touch "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_2.trim.fastq.gz"
              }
          )

          # Validate output files
          for f in "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_2.trim.fastq.gz"; do
            if [ ! -f "$f" ]; then
              echo "Creating empty file for failed sample $sample_name" >> trimming.log
              touch "$f"
            fi
            echo "Output file: $f ($(wc -c < "$f") bytes)" >> trimming.log
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
    continueOnReturnCode: true
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
      # Run FastQC on all files, continue even if some fail
      set +e
      fastqc -o qc_reports -t ~{cpu} ~{sep=' ' input_reads} 2>> qc.log
      set -e

      echo "Running MultiQC..." >> qc.log
      # Run MultiQC even if FastQC failed for some samples
      set +e
      multiqc qc_reports -o qc_reports --force 2>> qc.log
      set -e

      # Verify reports
      if [ $(ls qc_reports/*.{html,zip,txt} 2>/dev/null | wc -l) -eq 0 ]; then
        echo "WARNING: No QC reports generated, creating empty reports" >> qc.log
        touch qc_reports/empty_report.html
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
        sample_name=$(basename "$R1" | sed 's/_1.fastq.gz//; s/_1.trim.fastq.gz//')
        outdir="~{output_dir}/megahit_${sample_name}"

        echo "Assembling $sample_name (R1: $R1, R2: $R2)" >> assembly.log

        # Run assembly in a subshell to prevent failure from stopping the whole task
        (
          set +e
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
            contig_count=$(grep -c '^>' "~{output_dir}/${sample_name}_contigs.fa" || echo 0)
            echo "Generated contigs for $sample_name: $contig_count sequences" >> assembly.log
          else
            echo "WARNING: No contigs generated for $sample_name, creating empty file" >> assembly.log
            touch "~{output_dir}/${sample_name}_contigs.fa"
          fi
          set -e
        )
      done

      if [ $(find ~{output_dir} -name "*_contigs.fa" | wc -l) -eq 0 ]; then
        echo "ERROR: No contig files generated!" >> assembly.log
        exit 1
      fi

      echo "Assembly completed at $(date)" >> assembly.log
      echo "Output files created:" >> assembly.log
      find ~{output_dir} -type f -exec ls -lh {} \; >> assembly.log
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
          "$asm_file" || echo "PROKKA failed for $sample_name" >&2
      done
    else
      echo "Annotation skipped by user request" > annotation_results/skipped.txt
    fi
  >>>

  runtime {
    docker: "staphb/prokka:1.14.6"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
  }

  output {
    Array[File] annotation_output = if do_annotation then glob("annotation_results/*/*.gff") else []
    Array[String] annotation_dirs = if do_annotation then glob("annotation_results/*") else []
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

    mkdir -p gff_inputs
    for gff in ~{sep=' ' annotation_input}; do
      cp "$gff" gff_inputs/
    done

    echo "=== EXECUTING ROARY ===" >&2
    roary -f ~{output_prefix}_results \
          -p ~{cpu} \
          -e -n -v \
          -i 90 -cd 99 \
          gff_inputs/*.gff

    echo "=== LOCATING OUTPUT FILES ===" >&2
    output_dir="~{output_prefix}_results"
    core_alignment=$(find "$output_dir" -name "core_gene_alignment.aln" | head -1)

    if [ -z "$core_alignment" ]; then
      echo "ERROR: No core alignment file found" >&2
      exit 1
    fi

    mkdir -p final_output
    cp "$core_alignment" final_output/core_gene_alignment.aln
    cp "$output_dir"/gene_presence_absence.csv final_output/
    cp "$output_dir"/accessory_binary_genes.fa final_output/
    cp "$output_dir"/summary_statistics.txt final_output/
  >>>

  runtime {
    docker: "staphb/roary:3.13.0"
    memory: "8G"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: false
  }

  output {
    File gene_presence_absence = "final_output/gene_presence_absence.csv"
    File summary_statistics = "final_output/summary_statistics.txt"
    File accessory_binary = "final_output/accessory_binary_genes.fa"
    File core_alignment = "final_output/core_gene_alignment.aln"
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
    Int memory_gb = 16  # Increased from default
    Int disk_gb = 100   # Added explicit disk space
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

    # Check available memory
    echo "System memory information:" >> phylogeny.log
    free -h >> phylogeny.log

    if [ "~{do_phylogeny}" == "true" ]; then
      echo "Phylogeny parameters:" >> phylogeny.log
      echo "- Model: ~{model}" >> phylogeny.log
      echo "- CPU: ~{cpu}" >> phylogeny.log
      echo "- Memory allocated: ~{memory_gb} GB" >> phylogeny.log
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

      # Calculate expected memory needs based on alignment size
      alignment_size=$(wc -c < "~{alignment}")
      echo "Alignment size: $alignment_size bytes" >> phylogeny.log

      # Run FastTree with memory monitoring
      echo "Running FastTree for core genome..." >> phylogeny.log
      set +e
      ulimit -v $((~{memory_gb} * 1024 * 1024))  # Set memory limit

      FastTree ~{model} \
        -gamma \
        -quiet \
        -boot ~{bootstrap_replicates} \
        -log "phylogeny_results/core_~{tree_prefix}.log" \
        < "~{alignment}" > "phylogeny_results/core_~{tree_prefix}.nwk" 2>> phylogeny.log
      exit_code=$?
      set -e

      if [ $exit_code -ne 0 ]; then
        echo "WARNING: FastTree exited with code $exit_code" >> phylogeny.log

        # Attempt fallback with fewer bootstrap replicates if memory was the issue
        if grep -qi "oom" phylogeny.log || grep -qi "killed" phylogeny.log; then
          echo "Memory issue detected, retrying with 50 bootstrap replicates" >> phylogeny.log
          FastTree ~{model} \
            -gamma \
            -quiet \
            -boot 50 \
            -log "phylogeny_results/core_~{tree_prefix}_reduced.log" \
            < "~{alignment}" > "phylogeny_results/core_~{tree_prefix}.nwk" 2>> phylogeny.log || {
              echo "Fallback also failed, generating minimal tree" >> phylogeny.log
              echo "(A,B);" > "phylogeny_results/core_~{tree_prefix}.nwk"
            }
        else
          echo "Non-memory error, generating minimal tree" >> phylogeny.log
          echo "(A,B);" > "phylogeny_results/core_~{tree_prefix}.nwk"
        fi
      fi

      if [ ! -s "phylogeny_results/core_~{tree_prefix}.nwk" ]; then
        echo "ERROR: Core tree generation failed - empty output file" >> phylogeny.log
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
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2  # Allow preemptible instances
    continueOnReturnCode: true
    timeout: "24 hours"  # Extended timeout
  }

  output {
    File phylogeny_tree = if do_phylogeny then "phylogeny_results/core_~{tree_prefix}.nwk" else "phylogeny_results/core_~{tree_prefix}.nwk"
    File? phylogeny_log = if do_phylogeny then "phylogeny_results/core_~{tree_prefix}.log" else "skipped.txt"
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
      set +e
      FastTree ~{model} \
        -gamma \
        -quiet \
        -boot ~{bootstrap_replicates} \
        -log "phylogeny_results/accessory_~{tree_prefix}.log" \
        < "~{alignment}" > "phylogeny_results/accessory_~{tree_prefix}.nwk" 2>> phylogeny.log
      set -e

      if [ ! -s "phylogeny_results/accessory_~{tree_prefix}.nwk" ]; then
        echo "ERROR: Accessory tree generation failed - empty output file" >> phylogeny.log
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
    continueOnReturnCode: true
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

        # Run MLST in a subshell to prevent failure from stopping the whole task
        (
          set +e
          mlst "$asm_file" > "$output_file" 2>> mlst.log || {
            echo "MLST failed for $sample_name" >> mlst.log
            echo -e "file\tscheme\tst\trep1\trep2\trep3\trep4\trep5\trep6\trep7" > "$output_file"
          }
          set -e

          if [ -s "$output_file" ]; then
            awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
            mv "${output_file}.tmp" "$output_file"
          fi
        )
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

        (
          set +e
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
              --force 2>> ../variant.log || {
                echo "WARNING: Snippy failed for $SAMPLE_NAME" >> ../variant.log
              }
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
              --force 2>> ../variant.log || {
                echo "WARNING: Snippy failed for $SAMPLE_NAME" >> ../variant.log
              }
          fi
          set -e
        )

        cd ../..
      done

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
    File? local_db
    Boolean use_local_db = false
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
      echo "- Use local DB: ~{use_local_db}" >> amr.log
      echo "- Local DB path: ~{local_db}" >> amr.log
      echo "- Min identity: ~{minid}" >> amr.log
      echo "- Min coverage: ~{mincov}" >> amr.log
      echo "- CPU: ~{cpu}" >> amr.log

      mkdir -p amr_results

      if [ "~{use_local_db}" == "true" ] && [ -n "~{local_db}" ]; then
        echo "Setting up local AMR database" >> amr.log
        mkdir -p /root/abricate/db/resfinder_db
        cp "~{local_db}" /root/abricate/db/resfinder_db/resfinder.fa
        abricate --setupdb --db resfinder --debug >> amr.log 2>&1 || {
          echo "ERROR: Failed to setup local AMR database" >> amr.log
          exit 1
        }
        db_to_use="resfinder"
      else
        if ! abricate --list | grep -q "resfinder"; then
          abricate --setupdb --db resfinder >> amr.log 2>&1 || {
            echo "ERROR: Failed to setup resfinder database" >> amr.log
            exit 1
          }
        fi
        db_to_use="resfinder"
      fi

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
        output_file="amr_results/${sample_name}_amr.tsv"

        echo "Processing $sample_name" >> amr.log

        # Run AMR profiling in a subshell to prevent failure from stopping the whole task
        (
          set +e
          abricate \
            --db $db_to_use \
            --minid ~{minid} \
            --mincov ~{mincov} \
            --threads ~{cpu} \
            "$asm_file" > "$output_file" 2>> amr.log || {
              echo "AMR profiling failed for $sample_name" >> amr.log
              echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "$output_file"
            }
          set -e

          if [ -s "$output_file" ]; then
            awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
            mv "${output_file}.tmp" "$output_file"
          fi
        )
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
    File? local_db
    Boolean use_local_db = false
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
      echo "- Use local DB: ~{use_local_db}" >> plasmid.log
      echo "- Local DB path: ~{local_db}" >> plasmid.log
      echo "- CPU: ~{cpu}" >> plasmid.log

      mkdir -p plasmid_results

      if [ "~{use_local_db}" == "true" ] && [ -n "~{local_db}" ]; then
        echo "Setting up local plasmid database" >> plasmid.log
        mkdir -p /root/abricate/db/plasmidfinder
        cp "~{local_db}" /root/abricate/db/plasmidfinder/plasmidfinder.fa
        abricate --setupdb --db plasmidfinder --debug >> plasmid.log 2>&1 || {
          echo "ERROR: Failed to setup local plasmid database" >> plasmid.log
          exit 1
        }
        db_to_use="plasmidfinder"
      else
        if ! abricate --list | grep -q plasmidfinder; then
          abricate --setupdb --db plasmidfinder >> plasmid.log 2>&1 || {
            echo "ERROR: Failed to setup plasmidfinder database" >> plasmid.log
            exit 1
          }
        fi
        db_to_use="plasmidfinder"
      fi

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\..*//')
        echo "Processing $sample_name" >> plasmid.log

        # Run MGE analysis in a subshell to prevent failure from stopping the whole task
        (
          set +e
          abricate \
            --db $db_to_use \
            --mincov 80 \
            --minid 90 \
            --threads ~{cpu} \
            --nopath \
            "$asm_file" > "plasmid_results/${sample_name}.tsv" 2>> plasmid.log || {
              echo "Plasmid detection failed for $sample_name" >> plasmid.log
              echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "plasmid_results/${sample_name}.tsv"
            }
          set -e

          if [ -s "plasmid_results/${sample_name}.tsv" ]; then
            awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "plasmid_results/${sample_name}.tsv" > "plasmid_results/${sample_name}.tmp"
            mv "plasmid_results/${sample_name}.tmp" "plasmid_results/${sample_name}.tsv"
          fi
        )
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
    File? local_db
    Boolean use_local_db = false
    Int min_coverage = 60
    Float min_identity = 80.0
    Int cpu = 2
  }

  command <<<
    set -euo pipefail

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
    echo "- Use local DB: ~{use_local_db}" >> virulence.log
    echo "- Local DB path: ~{local_db}" >> virulence.log
    echo "- Min coverage: ~{min_coverage}" >> virulence.log
    echo "- Min identity: ~{min_identity}" >> virulence.log
    echo "- CPU: ~{cpu}" >> virulence.log

    if [ "~{use_local_db}" == "true" ] && [ -n "~{local_db}" ]; then
      echo "Setting up local virulence database" >> virulence.log
      mkdir -p /root/abricate/db/vfdb
      cp "~{local_db}" /root/abricate/db/vfdb/vfdb.fa
      abricate --setupdb --db vfdb --debug >> virulence.log 2>&1 || {
        echo "ERROR: Failed to setup local virulence database" >> virulence.log
        exit 1
      }
      db_to_use="vfdb"
    else
      if ! abricate --list | grep -q "~{db_name}"; then
        abricate --setupdb --db "~{db_name}" >> virulence.log 2>&1 || {
          echo "ERROR: Failed to setup ~{db_name} database" >> virulence.log
          exit 1
        }
      fi
      db_to_use="~{db_name}"
    fi

    mkdir -p virulence_results

    for asm_file in ~{sep=' ' assembly_output}; do
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//')
      output_file="virulence_results/${sample_name}_virulence.tsv"

      echo "Processing $sample_name" >> virulence.log

      (
        set +e
        abricate \
          --db $db_to_use \
          --mincov ~{min_coverage} \
          --minid ~{min_identity} \
          --nopath \
          --quiet \
          --threads ~{cpu} \
          "$asm_file" > "$output_file" 2>> virulence.log || {
            echo "Abricate failed for $sample_name" >> virulence.log
            echo -e "FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT" > "$output_file"
          }
        set -e

        if [ -s "$output_file" ]; then
          awk -v sample="$sample_name" 'BEGIN{OFS="\t"} NR==1 {print "SAMPLE",$0; next} {print sample,$0}' "$output_file" > "${output_file}.tmp"
          mv "${output_file}.tmp" "$output_file"
        fi
      )
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
    Array[File] virulence_reports = glob("virulence_results/*_virulence.tsv")
    File combined_report = "virulence_results/combined_virulence.tsv"
    File? virulence_log = "virulence.log"
  }
}

task BLAST_ANALYSIS {
  input {
    Array[File]+ contig_fastas
    String blast_db = "nt"
    File? local_blast_db
    Boolean use_local_blast = false
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

    function run_local_blast {
      local query="$1"
      local output_file="$2"
      local log_file="$3"

      echo "Attempting local BLAST..." >> "$log_file"

      # Check if we need to format the database
      if [ ! -f "~{local_blast_db}.nhr" ]; then
        echo "Formatting local BLAST database..." >> "$log_file"
        makeblastdb \
          -in "~{local_blast_db}" \
          -dbtype nucl \
          -out "~{local_blast_db}" \
          2>> "$log_file" || {
            echo "Failed to format local BLAST database" >> "$log_file"
            return 1
          }
      fi

      blastn \
        -query "$query" \
        -db "~{local_blast_db}" \
        -task blastn \
        -word_size 28 \
        -reward 1 -penalty -2 \
        -gapopen 2 -gapextend 1 \
        -outfmt "6 std qlen slen stitle" \
        -out "$output_file" \
        -evalue ~{evalue} \
        -max_target_seqs ~{max_target_seqs} \
        -num_threads ~{cpu} \
        2>> "$log_file"

      if [ -s "$output_file" ] && grep -q -v '^#' "$output_file"; then
        echo "Local BLAST succeeded" >> "$log_file"
        return 0
      else
        echo "Local BLAST failed or produced empty output" >> "$log_file"
        return 1
      fi
    }

    function run_remote_blast {
      local query="$1"
      local output_file="$2"
      local log_file="$3"
      local attempt=1
      local max_attempts=~{max_retries_per_sample}
      local delay=~{retry_delay_seconds}

      while [ $attempt -le $max_attempts ]; do
        echo "Attempt $attempt of $max_attempts for remote BLAST" >> "$log_file"

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
            echo "Remote BLAST succeeded on attempt $attempt" >> "$log_file"
            return 0
          else
            echo "Empty/invalid remote BLAST output on attempt $attempt" >> "$log_file"
          fi
        else
          echo "Remote BLAST failed with code $? on attempt $attempt" >> "$log_file"
        fi

        sleep $delay
        attempt=$((attempt + 1))
        delay=$((delay * 2))
      done

      echo "Max retries ($max_attempts) exceeded for remote BLAST" >> "$log_file"
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
      echo "- Local database: ~{local_blast_db}" >> blast.log
      echo "- Use local: ~{use_local_blast}" >> blast.log
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

        # First try local BLAST if configured
        local_success=false
        if [ "~{use_local_blast}" = true ] && [ -n "~{local_blast_db}" ]; then
          if run_local_blast "$filtered_contig" "$blast_output" "$blast_log"; then
            local_success=true
          fi
        fi

        # Fall back to remote BLAST if local failed or not configured
        if [ "$local_success" = false ]; then
          echo "Attempting remote BLAST as fallback" >> "$blast_log"
          run_remote_blast "$filtered_contig" "$blast_output" "$blast_log" || {
            echo "Both local and remote BLAST failed for $sample_id" >> "$blast_log"
            echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle" > "$blast_output"
          }
        fi

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
      echo "BLAST skipped by user request" > blast_results/skipped.txt
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

task TREE_VISUALIZATION {
  input {
    File input_tree
    Int width = 1200
    String image_format = "png"
    Int font_size = 8
    Boolean show_scale = true
    Boolean show_leaf_name = true
    Boolean show_support_values = true
    Int? continueOnReturnCode = 1
  }

  command <<<
    # Create all directories first
    mkdir -p /tmp/inputs /tmp/outputs final_phylogenetic_tree_image

    # Get the basename of the input file
    INPUT_BASENAME=$(basename "~{input_tree}")

    # Copy input file
    if [[ ! -f "~{input_tree}" ]]; then
      echo "Input file ~{input_tree} not found!" >&2
      exit 1
    fi
    cp "~{input_tree}" /tmp/inputs/input.nwk

    # Execute Python script
    python3 <<'PYTHON_SCRIPT' > /tmp/outputs/render.log 2>&1
import os
import sys
import traceback
from ete3 import Tree, TreeStyle, TextFace, NodeStyle

def main():
    try:
        print("[DEBUG] Starting tree rendering")
        input_path = "/tmp/inputs/input.nwk"
        output_path = f"/tmp/outputs/tree.~{image_format}"

        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file missing: {input_path}")

        t = Tree(input_path, format=1)
        print(f"[SUCCESS] Loaded tree with {len(t)} leaves")

        ts = TreeStyle()
        ts.show_scale = ~{true="True" false="False" show_scale}
        ts.scale_length = 0.1
        ts.mode = "r"
        ts.rotation = 0
        ts.branch_vertical_margin = 15
        ts.show_leaf_name = ~{true="True" false="False" show_leaf_name}
        ts.min_leaf_separation = 5
        ts.allow_face_overlap = True
        ts.complete_branch_lines_when_necessary = True
        ts.root_opening_factor = 0.5
        ts.margin_left = 50
        ts.margin_right = 50
        ts.margin_top = 50
        ts.margin_bottom = 50

        # Add leaf names with controlled font size if not showing by default
        if not ~{true="True" false="False" show_leaf_name}:
            for leaf in t.iter_leaves():
                face = TextFace(leaf.name, fsize=~{font_size})
                leaf.add_face(face, column=0, position="branch-right")

        # Add support values if enabled
        if ~{true="True" false="False" show_support_values}:
            for node in t.traverse():
                if not node.is_leaf():
                    ns = NodeStyle()
                    ns["size"] = 0
                    if hasattr(node, 'support'):
                        support_value = node.support
                        if support_value is not None:
                            support_text = f"{support_value:.0%}"
                            support_face = TextFace(support_text,
                                                 fsize=~{font_size},
                                                 fgcolor="black",
                                                 bold=True)
                            support_face.margin_right = 10
                            support_face.margin_left = 10
                            node.add_face(support_face, column=0, position="branch-top")
                    node.set_style(ns)

        print(f"[DEBUG] Rendering to {output_path}")
        t.render(output_path, w=~{width}, units="px", tree_style=ts)

        if not os.path.exists(output_path):
            raise RuntimeError("Rendering completed but no output file created")

        print("[SUCCESS] Render completed")
        return 0

    except Exception as e:
        print(f"[ERROR] {traceback.format_exc()}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
PYTHON_SCRIPT

    # Handle outputs with new directory structure
    if [[ -f "/tmp/outputs/tree.~{image_format}" ]]; then
      mkdir -p final_phylogenetic_tree_image
      cp "/tmp/outputs/tree.~{image_format}" "final_phylogenetic_tree_image/phylogenetic_tree_${INPUT_BASENAME}.~{image_format}"
      cp "/tmp/outputs/render.log" "./render_${INPUT_BASENAME}.log"
      exit 0
    else
      mkdir -p final_phylogenetic_tree_image
      echo "Rendering failed" > "error_${INPUT_BASENAME}.log"
      cp "/tmp/outputs/render.log" "./render_${INPUT_BASENAME}.log"
      exit 1
    fi
  >>>

  runtime {
    docker: "gmboowa/ete3-render:1.14"
    memory: "4 GB"
    cpu: 2
    continueOnReturnCode: true
  }

  output {
    File? final_image = "final_phylogenetic_tree_image/phylogenetic_tree_~{basename(input_tree)}.~{image_format}"
    File render_log = "render_~{basename(input_tree)}.log"
  }
}
