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
    Int max_cpus = 8
    Int max_memory_gb = 16
    Int min_assembly_quality = 50
    Int min_read_length = 50
    Int min_mapping_quality = 20
    Int tree_image_width = 1200
    String tree_image_format = "png"
    Int tree_font_size = 8
    File? logo_file
    Boolean do_virulence = true
    # Derived values with adjusted resources
    Int cpu_4 = if (max_cpus < 4) then max_cpus else 4
    Int cpu_8 = if (max_cpus < 8) then max_cpus else 8
    Int cpu_16 = max_cpus
    Int cpu_2 = if (max_cpus < 2) then max_cpus else 2
    Int mem_8 = if (max_memory_gb < 8) then max_memory_gb else 8
    Int mem_12 = if (max_memory_gb < 12) then max_memory_gb else 12
    Int mem_16 = max_memory_gb
  }

  meta {
    workflow_timeout: "168 hours"
    workflow_heartbeat_interval: "10 minutes"
    workflow_heartbeat_ttl: "30 minutes"
    allowNestedInputs: true
    maxRetries: 3
    continueOnReturnCode: "0,1"
    author: "Gerald Mboowa, Ivan Sserwadda & Stephen Kanyerezi"
    email: "gmboowa@gmail.com | ivangunz23@gmail.com | kanyerezi30@gmail.com"
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
      cpu = cpu_8,
      memory_gb = mem_16
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
      cpu = cpu_8,
      memory_gb = mem_12,
      min_quality = min_mapping_quality
  }

  call PANGENOME {
    input:
      annotation_input = ANNOTATION.annotation_output,
      do_pangenome = do_pangenome,
      cpu = cpu_8,
      memory_gb = mem_16,
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

  # >>> REPLACE subworkflow call with a direct scatter over the MGE task <<<
  scatter (asm in ASSEMBLY.assembly_output) {
    call MGE {
      input:
        assembly        = asm,
        do_mge_analysis = do_mge_analysis,
        local_db        = local_mge_db,
        use_local_db    = use_local_blast,
        cpu             = cpu_4
    }
  }

  if (do_virulence) {
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
      cpu = cpu_8,
      memory_gb = mem_12
  }

  if (do_phylogeny) {
    call CORE_PHYLOGENY {
      input:
        alignment = PANGENOME.core_alignment,
        do_phylogeny = true,
        model = phylogeny_model,
        cpu = cpu_8,
        memory_gb = mem_16,
        tree_prefix = "core_genes",
        bootstrap_replicates = 100
    }

    call ACCESSORY_PHYLOGENY {
      input:
        alignment = PANGENOME.accessory_binary,
        do_phylogeny = true,
        model = phylogeny_model,
        cpu = cpu_8,
        memory_gb = mem_16,
        tree_prefix = "accessory_genes",
        bootstrap_replicates = 100
    }

    if (size(CORE_PHYLOGENY.phylogeny_tree, "B") > 10.0) {
      call TREE_VISUALIZATION as CORE_TREE {
        input:
          input_tree   = CORE_PHYLOGENY.phylogeny_tree,
          width        = tree_image_width,
          image_format = tree_image_format,
          font_size    = tree_font_size,
          tree_title   = "Core Genes Phylogenetic Tree"
      }
    }

    if (defined(ACCESSORY_PHYLOGENY.phylogeny_tree) && size(ACCESSORY_PHYLOGENY.phylogeny_tree, "B") > 10.0) {
      call TREE_VISUALIZATION as ACCESSORY_TREE {
        input:
          input_tree   = ACCESSORY_PHYLOGENY.phylogeny_tree,
          width        = tree_image_width,
          image_format = tree_image_format,
          font_size    = tree_font_size,
          tree_title   = "Accessory Genes Phylogenetic Tree"
      }
    }
  }

  Array[File] empty_files = []

  if (do_reporting) {
    call MERGE_REPORTS {
      input:
        # Core sections
        quality_reports          = QUALITY_CONTROL.quality_reports,
        trimming_report_html     = TRIMMING.trimming_report,
        assembly_stats_html      = ASSEMBLY.assembly_stats,
        annotation_summary_html  = ANNOTATION.annotation_summary,
        pangenome_report_html    = PANGENOME.pangenome_report,
        gene_heatmap_png         = PANGENOME.gene_heatmap,
        pangenome_report         = PANGENOME.pangenome_report,
        gene_heatmap             = PANGENOME.gene_heatmap,
        mlst_combined_html       = MLST.combined_html,
        variant_summary_html     = VARIANT_CALLING.variant_summary,

        # Per-sample sections
        amr_reports              = AMR_PROFILING.html_reports,
        mge_reports              = MGE.html_out,
        virulence_reports        = if do_virulence then VIRULENCE_ANALYSIS.html_reports else empty_files,

        # BLAST per-sample HTMLs
        blast_reports            = BLAST_ANALYSIS.blast_reports,

        # Trees (may be empty depending on upstream conditionals)
        tree_images              = select_all([CORE_TREE.final_image, ACCESSORY_TREE.final_image]),

        # Meta
        workflow_name            = "rMAP Analysis Pipeline",
        version                  = "2.0",
        footer_sentence          = "This report was generated by rMAP 2.0 pipeline",
        logo_file                = logo_file,
        skip                     = false
    }
  }

  output {
    Array[File] quality_reports = QUALITY_CONTROL.quality_reports
    Array[File] trimmed_reads   = TRIMMING.trimmed_reads
    Array[File] assembly_output = ASSEMBLY.assembly_output
    String      assembly_dir    = ASSEMBLY.assembly_dir_out
    Array[File]   annotation_output = ANNOTATION.annotation_output
    Array[String] annotation_dirs   = ANNOTATION.annotation_dirs
    Array[File] mlst_output   = MLST.mlst_outputs
    File        combined_mlst = MLST.combined_mlst
    Array[File]   variant_output = VARIANT_CALLING.variant_output
    Array[File]   variant_dirs   = VARIANT_CALLING.variant_dirs
    String        variants_dir   = VARIANT_CALLING.variants_dir
    File gene_presence_absence = PANGENOME.gene_presence_absence
    File core_alignment        = PANGENOME.core_alignment
    File accessory_alignment   = PANGENOME.accessory_binary
    File? core_phylogeny_output      = CORE_PHYLOGENY.phylogeny_tree
    File? accessory_phylogeny_output = ACCESSORY_PHYLOGENY.phylogeny_tree
    File? core_tree_image           = CORE_TREE.final_image
    File? accessory_tree_image      = ACCESSORY_TREE.final_image
    File? core_tree_render_log      = CORE_TREE.render_log
    File? accessory_tree_render_log = ACCESSORY_TREE.render_log
    Array[File] amr_output   = AMR_PROFILING.amr_outputs
    File        combined_amr = AMR_PROFILING.combined_amr
    Array[File] plasmid_results       = MGE.tsv_out
    Array[File] plasmid_html_reports  = MGE.html_out
    Array[File] blast_results = BLAST_ANALYSIS.blast_results
    Array[File] blast_top10   = BLAST_ANALYSIS.blast_top10
    Array[File] blast_logs    = BLAST_ANALYSIS.blast_logs
    Array[File]? virulence_reports         = VIRULENCE_ANALYSIS.virulence_reports
    File? final_report_html = MERGE_REPORTS.final_report_html
    File? final_report_tgz  = MERGE_REPORTS.final_report_tgz
  }
}


task CONFIGURATION {
  input {
    Int max_cpus = 16
    Int max_memory_gb = 16
  }

  command <<<
    echo "System Configuration:" > config.log
    echo "Max CPUs: ~{max_cpus}" >> config.log
    echo "Max Memory: ~{max_memory_gb} GB" >> config.log
    echo "Workflow version: 2.0" >> config.log
    echo "Start time: $(date)" >> config.log
    echo "Hostname: $(hostname)" >> config.log
    echo "Kernel: $(uname -a)" >> config.log
    cat config.log
  >>>

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
    Boolean skip = false
  }

  command <<<
    set -euo pipefail

    # Function to calculate read statistics
    calculate_read_stats() {
      local file=$1
      local stats_file=$2
      local read_end=$3

      if [ ! -s "$file" ]; then
        echo "Metric,${read_end}" >> "$stats_file"
        echo "Number of reads,0" >> "$stats_file"
        echo "Average length,0" >> "$stats_file"
        echo "GC content,0" >> "$stats_file"
        return
      fi

      local num_reads=$(gzip -cd -- "$file" | awk 'END {print NR/4}')
      local avg_length=$(gzip -cd -- "$file" | awk 'NR%4==2 {sum+=length($0)} END {print sum/(NR/4)}')
      local gc_content=$(gzip -cd -- "$file" | awk 'NR%4==2 {gc+=gsub(/[GC]/, ""); at+=gsub(/[AT]/, "")} END {print (gc/(gc+at))*100}')

      echo "Metric,${read_end}" >> "$stats_file"
      echo "Number of reads,$num_reads" >> "$stats_file"
      echo "Average length,$avg_length" >> "$stats_file"
      echo "GC content,$gc_content" >> "$stats_file"
    }

    if [ "~{skip}" == "true" ]; then
      echo "Skipping trimming process as requested" > trimming_skipped.log
      echo "Creating empty output files for portability" >> trimming_skipped.log

      mkdir -p trimmed
      mkdir -p read_stats

      counter=0
      R1=""
      for file in ~{sep=' ' input_reads}; do
        if [ $((counter % 2)) -eq 0 ]; then
          R1="$file"
        else
          R2="$file"
          R1_base=$(basename "$R1")
          sample_name=$(echo "$R1_base" | sed -e 's/[._][Rr]1[._].*//' -e 's/[._]1[._].*//')

          touch "trimmed/${sample_name}_1.trim.fastq.gz"
          touch "trimmed/${sample_name}_2.trim.fastq.gz"
          touch "read_stats/${sample_name}_stats.csv"
          echo "Created empty output files for sample: $sample_name" >> trimming_skipped.log
        fi
        counter=$((counter + 1))
      done

      echo "Sample,Input_Pairs,Both_Surviving,Forward_Only,Reverse_Only,Dropped" > trimming_stats.csv
      echo "0,0,0,0,0" >> trimming_stats.csv

      cat > trimming_report.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>Trimming Statistics Report - Skipped</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .skip-notice {
            background-color: #fff3cd;
            padding: 20px;
            border-radius: 5px;
            margin: 20px 0;
            text-align: center;
        }
    </style>
</head>
<body>
    <div class="skip-notice">
        <h2>Trimming Process Skipped</h2>
        <p>Read trimming was skipped as requested in the workflow parameters.</p>
        <p>Generated: $(date)</p>
    </div>
</body>
</html>
EOF

      echo "Skipping completed at $(date)" >> trimming_skipped.log
      exit 0
    fi

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
      mkdir -p read_stats
      echo "Trimming parameters:" >> trimming.log
      echo "- CPU: ~{cpu}" >> trimming.log
      echo "- Min length: ~{min_length}" >> trimming.log

      echo "Sample,Input_Pairs,Both_Surviving,Forward_Only,Reverse_Only,Dropped" > trimming_stats.csv

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

          sample_log="trimmed/${sample_name}.log"
          touch "$sample_log"

          (
            trimmomatic PE -threads ~{cpu} \
              "$R1" "$R2" \
              "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_1.unpair.fastq.gz" \
              "trimmed/${sample_name}_2.trim.fastq.gz" "trimmed/${sample_name}_2.unpair.fastq.gz" \
              ILLUMINACLIP:~{adapters}:2:30:10:8:true \
              LEADING:20 \
              TRAILING:20 \
              SLIDINGWINDOW:4:20 \
              MINLEN:~{min_length} > "$sample_log" 2>&1 || {
                echo "WARNING: Trimmomatic failed for sample $sample_name" >> trimming.log
                touch "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_2.trim.fastq.gz"
              }
          )

          stats_line=$(grep "Input Read Pairs" "$sample_log" || true)
          if [[ -n "${stats_line:-}" ]]; then
            input_pairs=$(awk '{print $4}' <<< "$stats_line")
            both=$(awk '{print $7}' <<< "$stats_line")
            forward=$(awk '{print $12}' <<< "$stats_line")
            reverse=$(awk '{print $17}' <<< "$stats_line")
            dropped=$(awk '{print $20}' <<< "$stats_line")

            echo "$sample_name,$input_pairs,$both,$forward,$reverse,$dropped" >> trimming_stats.csv
          else
            echo "$sample_name,0,0,0,0,0" >> trimming_stats.csv
            echo "WARNING: Could not extract statistics for sample $sample_name" >> trimming.log
          fi

          > "read_stats/${sample_name}_stats.csv"

          if [ -s "trimmed/${sample_name}_1.trim.fastq.gz" ]; then
            calculate_read_stats "trimmed/${sample_name}_1.trim.fastq.gz" \
              "read_stats/${sample_name}_stats.csv" "R1"
          fi

          if [ -s "trimmed/${sample_name}_2.trim.fastq.gz" ]; then
            calculate_read_stats "trimmed/${sample_name}_2.trim.fastq.gz" \
              "read_stats/${sample_name}_stats.csv" "R2"
          fi

          cat "$sample_log" >> trimming.log
          rm "$sample_log"

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

      cat > trimming_report.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>Trimming Statistics Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }
        .good { color: #27ae60; }
        .warning { color: #f39c12; }
        .bad { color: #e74c3c; }
        .stats-table { margin-top: 30px; }
    </style>
</head>
<body>
    <h1>Trimming Statistics Report</h1>
    <div class="summary">
        <p><strong>Date:</strong> $(date)</p>
        <p><strong>Parameters:</strong></p>
        <ul>
            <li>Minimum length: ~{min_length} bp</li>
            <li>Threads: ~{cpu}</li>
        </ul>
    </div>

    <h2>Trimming Statistics</h2>
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Input Pairs</th>
                <th>Both Surviving</th>
                <th>Forward Only</th>
                <th>Reverse Only</th>
                <th>Dropped (%)</th>
                <th>Survival Rate (%)</th>
            </tr>
        </thead>
        <tbody>
EOF

      while IFS=',' read -r sample input both forward reverse dropped; do
        if [ "$sample" == "Sample" ]; then
          continue
        fi

        if [ "$input" -gt 0 ]; then
          total_surviving=$((both + forward + reverse))
          survival_rate=$(awk -v ts="$total_surviving" -v i="$input" 'BEGIN {printf "%.2f", ts*100/i}')
          dropped_rate=$(awk -v d="$dropped" -v i="$input" 'BEGIN {printf "%.2f", d*100/i}')
        else
          survival_rate="0.00"
          dropped_rate="0.00"
        fi

        if (( $(awk -v sr="$survival_rate" 'BEGIN {print (sr >= 90)}') )); then
          class="good"
        elif (( $(awk -v sr="$survival_rate" 'BEGIN {print (sr >= 70)}') )); then
          class="warning"
        else
          class="bad"
        fi

        echo "            <tr>" >> trimming_report.html
        echo "                <td>$sample</td>" >> trimming_report.html
        echo "                <td>$input</td>" >> trimming_report.html
        echo "                <td>$both</td>" >> trimming_report.html
        echo "                <td>$forward</td>" >> trimming_report.html
        echo "                <td>$reverse</td>" >> trimming_report.html
        echo "                <td>${dropped_rate}%</td>" >> trimming_report.html
        echo "                <td class=\"$class\">${survival_rate}%</td>" >> trimming_report.html
        echo "            </tr>" >> trimming_report.html
      done < trimming_stats.csv

      cat >> trimming_report.html <<EOF
        </tbody>
    </table>
</body>
</html>
EOF

      echo "Trimming completed successfully at $(date)" >> trimming.log
      echo "Output files created:" >> trimming.log
      ls -lh trimmed/* >> trimming.log
      ls -lh read_stats/* >> trimming.log

      echo "Trimming executed (not skipped) on $(date)" > trimming_skipped.log
    else
      echo "Trimming skipped by user request" > trimming_skipped.txt
      echo "Trimming skipped on $(date)" > trimming_skipped.log
    fi
  >>>

  runtime {
    docker: "staphb/trimmomatic:0.39"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: [0, -1]
    preemptible: 2
    timeout: "8 hours"
  }

  output {
    Array[File] trimmed_reads = if do_trimming then glob("trimmed/*_[12].trim.fastq.gz") else []
    Array[File]? read_stats = if do_trimming then glob("read_stats/*_stats.csv") else []
    File? trimming_log = if do_trimming then "trimming.log" else "trimming_skipped.txt"
    File? trimming_stats = if do_trimming then "trimming_stats.csv" else "trimming_skipped.txt"
    File? trimming_report = if do_trimming then "trimming_report.html" else "trimming_skipped.txt"
    File? skip_log = "trimming_skipped.log"
  }
}
task QUALITY_CONTROL {
  input {
    Array[File]+ input_reads
    Boolean do_quality_control
    Int cpu = 4
    Boolean run_multiqc = true
    Boolean skip = false
  }

  command <<<
    set -euo pipefail

    if [ "~{skip}" == "true" ]; then
      echo "Skipping quality control as requested" > qc_skip.log
      echo "Creating minimal output files for portability" >> qc_skip.log

      mkdir -p qc_reports

      # Create empty output files
      touch qc_reports/summary.html
      touch qc_reports/skipped.txt

      # Create minimal HTML report
      cat > qc_reports/summary.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>Quality Control Report - Skipped</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .skip-notice {
            background-color: #fff3cd;
            padding: 20px;
            border-radius: 5px;
            margin: 20px 0;
            text-align: center;
        }
    </style>
</head>
<body>
    <div class="skip-notice">
        <h2>Quality Control Skipped</h2>
        <p>Quality control analysis was skipped as requested in the workflow parameters.</p>
        <p>Generated: $(date)</p>
    </div>
</body>
</html>
EOF

      echo "Skipping completed at $(date)" >> qc_skip.log
      exit 0
    fi

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
      echo "- Run MultiQC: ~{run_multiqc}" >> qc.log

      echo "Running FastQC..." >> qc.log
      # Run FastQC on all files, continue even if some fail
      set +e
      fastqc -o qc_reports -t ~{cpu} ~{sep=' ' input_reads} 2>> qc.log
      set -e

      # Generate individual report links file with consistent naming
      cat > qc_reports/individual_reports.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>Individual QC Reports</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h2>Individual QC Reports</h2>
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Report</th>
                <th>Raw Data</th>
            </tr>
        </thead>
        <tbody>
EOF

      # Add entries for each FastQC report with consistent naming
      idx=1
      for f in qc_reports/*_fastqc.html; do
        sample=$(basename "$f" | sed 's/_fastqc\.html//')
        report_name="QC_Report_${idx}_${sample}.html"
        cp "$f" "qc_reports/${report_name}"
        echo "            <tr>" >> qc_reports/individual_reports.html
        echo "                <td>${sample}</td>" >> qc_reports/individual_reports.html
        echo "                <td><a href=\"${report_name}\">QC Report ${idx} — ${sample}</a></td>" >> qc_reports/individual_reports.html
        zip_file="${f%.html}.zip"
        new_zip_name="QC_Data_${idx}_${sample}.zip"
        cp "$zip_file" "qc_reports/${new_zip_name}"
        echo "                <td><a href=\"${new_zip_name}\">Download Data</a></td>" >> qc_reports/individual_reports.html
        echo "            </tr>" >> qc_reports/individual_reports.html
        idx=$((idx+1))
      done

      cat >> qc_reports/individual_reports.html <<EOF
        </tbody>
    </table>
</body>
</html>
EOF

      if [ "~{run_multiqc}" == "true" ]; then
        echo "Running MultiQC..." >> qc.log

        # Create a simple MultiQC report alternative if MultiQC isn't available
        if ! command -v multiqc &> /dev/null; then
          echo "WARNING: MultiQC not found, generating simplified report" >> qc.log

          # Generate a simple HTML summary of FastQC results
          cat > qc_reports/fastqc_summary.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>QC Summary Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>QC Summary Report</h1>
    <p>Generated at $(date)</p>
    <p><a href="individual_reports.html">View Individual Reports</a></p>
EOF

          # Try to extract some basic metrics from FastQC data
          idx=1
          for f in qc_reports/*_fastqc/fastqc_data.txt; do
            if [ -f "$f" ]; then
              sample=$(basename "$f" | sed 's/_fastqc\.fastqc_data\.txt//')
              basic_stats=$(awk '/>>Basic Statistics/,/>>END_MODULE/' "$f" | grep -v '>>')

              echo "<h3>QC Report ${idx} — ${sample}</h3>" >> qc_reports/fastqc_summary.html
              echo "<table>" >> qc_reports/fastqc_summary.html
              echo "$basic_stats" | while read -r line; do
                if [[ "$line" != "" ]]; then
                  key=$(echo "$line" | cut -f1)
                  value=$(echo "$line" | cut -f2)
                  echo "<tr><td>$key</td><td>$value</td></tr>" >> qc_reports/fastqc_summary.html
                fi
              done
              echo "</table>" >> qc_reports/fastqc_summary.html
              idx=$((idx+1))
            fi
          done

          cat >> qc_reports/fastqc_summary.html <<EOF
    <p>Note: Full MultiQC report not available. Install MultiQC for more comprehensive analysis.</p>
</body>
</html>
EOF
        else
          # Run MultiQC if available
          multiqc qc_reports -o qc_reports --force 2>> qc.log || {
            echo "ERROR: MultiQC failed" >> qc.log
            # Create fallback summary
            cp qc_reports/individual_reports.html qc_reports/fastqc_summary.html
          }
        fi
      fi

      # Verify reports
      if [ $(ls qc_reports/*.{html,zip} 2>/dev/null | wc -l) -eq 0 ]; then
        echo "WARNING: No QC reports generated, creating empty reports" >> qc.log
        touch qc_reports/empty_report.html
      fi

      # Create a canonical QC summary file we can output consistently
      if [ -f qc_reports/multiqc_report.html ]; then
        cp qc_reports/multiqc_report.html qc_reports/summary.html
      elif [ -f qc_reports/fastqc_summary.html ]; then
        cp qc_reports/fastqc_summary.html qc_reports/summary.html
      elif [ -f qc_reports/individual_reports.html ]; then
        cp qc_reports/individual_reports.html qc_reports/summary.html
      elif [ -f qc_reports/empty_report.html ]; then
        cp qc_reports/empty_report.html qc_reports/summary.html
      else
        touch qc_reports/summary.html
      fi

      echo "Quality control completed at $(date)" >> qc.log
      echo "Output files created:" >> qc.log
      ls -lh qc_reports/* >> qc.log

      # Ensure skip marker exists even when not skipping, so outputs never break
      echo "QC executed (not skipped) on $(date)" > qc_skip.log
    else
      mkdir -p qc_reports
      echo "QC skipped by user request" > qc_reports/skipped.txt
      # Also ensure skip marker exists for outputs
      echo "QC skipped on $(date)" > qc_skip.log
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
    Array[File] quality_reports = if do_quality_control then glob("qc_reports/QC_Report_*.html") else ["qc_reports/skipped.txt"]
    File individual_reports = if do_quality_control then "qc_reports/individual_reports.html" else "qc_reports/skipped.txt"
    File? qc_log = if do_quality_control then "qc.log" else "qc_reports/skipped.txt"
    File fastqc_summary_html = if do_quality_control then "qc_reports/summary.html" else "qc_reports/skipped.txt"
    File skip_log = "qc_skip.log"
  }
}


task ASSEMBLY {
  input {
    Array[File]+ input_reads
    String assembler = "megahit"
    Boolean do_assembly = true
    String output_dir = "assembly"
    Int cpu = 8
    Int memory_gb = 16
    Int min_quality = 50
    Boolean run_quast = false
    Float memory_usage = 0.8
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
      echo "- Memory usage fraction: ~{memory_usage}" >> assembly.log
      echo "- Min quality: ~{min_quality}" >> assembly.log
      echo "- Run QUAST: ~{run_quast}" >> assembly.log

      files=( ~{sep=' ' input_reads} )
      if [ $(( ${#files[@]} % 2 )) -ne 0 ]; then
          echo "ERROR: Odd number of input files for assembly" >> assembly.log
          exit 1
      fi

      # Create HTML report header
      cat > ~{output_dir}/assembly_stats.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>Assembly Statistics Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>Assembly Statistics Report</h1>
    <p>Generated at $(date)</p>
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Contigs</th>
                <th>Total Length</th>
                <th>Min Length</th>
                <th>Max Length</th>
                <th>Avg Length</th>
                <th>N50</th>
            </tr>
        </thead>
        <tbody>
EOF

      for ((i=0; i<"${#files[@]}"; i+=2)); do
        R1="${files[i]}"
        R2="${files[i+1]}"
        sample_name=$(basename "$R1" | sed 's/_1.fastq.gz//; s/_1.trim.fastq.gz//; s/_R1.*//; s/_1.*//')
        outdir="~{output_dir}/megahit_${sample_name}"
        contig_file="~{output_dir}/${sample_name}.fa"  # Removed _contigs suffix

        echo "Assembling $sample_name (R1: $R1, R2: $R2)" >> assembly.log

        # Run assembly in a subshell
        (
          set +e
          megahit \
            -1 "$R1" -2 "$R2" \
            -o "$outdir" \
            -t ~{cpu} \
            --memory ~{memory_usage} \
            --min-count 2 \
            --min-contig-len ~{min_quality} \
            --k-list 21,33,55,77 \
            --merge-level 20,0.98 \
            --prune-level 2 \
            --prune-depth 2 \
            --no-mercy \
            2>> assembly.log

          if [ -f "$outdir/final.contigs.fa" ]; then
            cp "$outdir/final.contigs.fa" "$contig_file"

            # Extract statistics from MEGAHIT log
            stats=$(grep -A 1 "Merging to output final contigs" "$outdir/log" | tail -1 | \
                   sed -n 's/.* \([0-9]\+\) contigs, total \([0-9]\+\) bp, min \([0-9]\+\) bp, max \([0-9]\+\) bp, avg \([0-9]\+\) bp, N50 \([0-9]\+\) bp.*/\1 \2 \3 \4 \5 \6/p')

            if [ -n "$stats" ]; then
              contigs=$(echo $stats | awk '{print $1}')
              total_len=$(echo $stats | awk '{print $2}')
              min_len=$(echo $stats | awk '{print $3}')
              max_len=$(echo $stats | awk '{print $4}')
              avg_len=$(echo $stats | awk '{print $5}')
              n50=$(echo $stats | awk '{print $6}')

              echo "Generated contigs for $sample_name: $contigs sequences" >> assembly.log

              # Add to HTML report
              echo "<tr><td>$sample_name</td><td>$contigs</td><td>$total_len</td><td>$min_len</td><td>$max_len</td><td>$avg_len</td><td>$n50</td></tr>" >> ~{output_dir}/assembly_stats.html
            else
              echo "WARNING: Could not extract stats for $sample_name from log" >> assembly.log
              echo "<tr><td>$sample_name</td><td colspan='6'>Statistics unavailable</td></tr>" >> ~{output_dir}/assembly_stats.html
            fi
          else
            echo "WARNING: No contigs generated for $sample_name, creating empty file" >> assembly.log
            touch "$contig_file"
            echo "<tr><td>$sample_name</td><td colspan='6'>Assembly failed</td></tr>" >> ~{output_dir}/assembly_stats.html
          fi
          set -e
        )
      done

      # Close HTML report
      cat >> ~{output_dir}/assembly_stats.html <<EOF
        </tbody>
    </table>
    <p>Note: For more detailed analysis, consider running QUAST separately.</p>
</body>
</html>
EOF

      if [ $(find ~{output_dir} -name "*.fa" | wc -l) -eq 0 ]; then
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
    continueOnReturnCode: [0, -9]
    timeout: "24 hours"
  }

  output {
    Array[File] assembly_output = if do_assembly then glob("~{output_dir}/*.fa") else []  # Changed pattern
    File assembly_stats = if do_assembly then "~{output_dir}/assembly_stats.html" else "~{output_dir}/skipped.txt"
    String assembly_dir_out = "~{output_dir}"
    File? assembly_log = if do_assembly then "assembly.log" else "~{output_dir}/skipped.txt"
  }
}

task ANNOTATION {
  input {
    Array[File]+ assembly_output
    Boolean do_annotation
    Int cpu = 8
    Int memory_gb = 16
    Boolean skip = false
  }

  command <<<
    # Always create a status/skip log so the output is never missing
    : > skip_annotation.log
    echo "ANNOTATION started on $(date)" >> skip_annotation.log

    if [ "~{skip}" == "true" ]; then
      echo "Skipping annotation as requested" >> skip_annotation.log
      mkdir -p annotation_results
      touch annotation_results/skipped.txt

      cat > annotation_results/summary.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>Annotation Summary - Skipped</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .skip-notice {
            background-color: #fff3cd;
            padding: 20px;
            border-radius: 5px;
            margin: 20px 0;
            text-align: center;
        }
    </style>
</head>
<body>
    <div class="skip-notice">
        <h2>Annotation Process Skipped</h2>
        <p>Genome annotation was skipped as requested in the workflow parameters.</p>
        <p>Generated: $(date)</p>
    </div>
</body>
</html>
EOF
      echo "Skipping completed at $(date)" >> skip_annotation.log
      exit 0
    fi

    if [ "~{do_annotation}" == "true" ]; then
      mkdir -p annotation_results

      # Initialize HTML report
      echo "<html>
      <head>
        <title>Annotation Summary</title>
        <style>
          body { font-family: Arial, sans-serif; margin: 20px; }
          table { border-collapse: collapse; width: 50%; margin-top: 20px; }
          th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
          th { background-color: #3498db; color: white; }
          tr:nth-child(even) { background-color: #f2f2f2; }
        </style>
      </head>
      <body>
        <h1>Annotation Summary Report</h1>
        <table>
          <tr><th>Sample ID</th><th>Gene Count</th></tr>" > annotation_results/summary.html

      for asm_file in ~{sep=' ' assembly_output}; do
        sample_name=$(basename "$asm_file" | sed 's/\..*//')
        output_dir="annotation_results/${sample_name}"

        echo "Running PROKKA annotation for sample: $sample_name" >> skip_annotation.log

        # Run PROKKA
        prokka \
          --outdir "$output_dir" \
          --prefix "$sample_name" \
          --cpus ~{cpu} \
          --mincontiglen 200 \
          --compliant \
          --force \
          "$asm_file" || {
            echo "PROKKA failed for $sample_name" >&2
            # Add failed sample to report
            echo "<tr><td>$sample_name</td><td>FAILED</td></tr>" >> annotation_results/summary.html
            continue
          }

        # Count genes in the GFF file
        gene_count=0
        if [ -f "$output_dir/$sample_name.gff" ]; then
          gene_count=$(grep -c $'\tCDS\t' "$output_dir/$sample_name.gff" || echo 0)
        fi

        # Add to HTML report
        echo "<tr><td>$sample_name</td><td>$gene_count</td></tr>" >> annotation_results/summary.html
      done

      # Close HTML report
      echo "</table>
        <p>Report generated on $(date)</p>
      </body>
      </html>" >> annotation_results/summary.html

    else
      echo "Annotation skipped by user request" > annotation_results/skipped.txt
    fi
  >>>

  runtime {
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 50 HDD"
    docker: "staphb/prokka:1.14.5"
    continueOnReturnCode: true
  }

  output {
    Array[File] annotation_output = if do_annotation then glob("annotation_results/*/*.gff") else []
    Array[String] annotation_dirs = if do_annotation then glob("annotation_results/*") else []
    File? annotation_summary = if do_annotation then "annotation_results/summary.html" else "annotation_results/skipped.txt"
    File skip_log = "skip_annotation.log"
  }
}

task PANGENOME {
  input {
    Array[File] annotation_input
    Boolean do_pangenome = true
    String output_prefix = "pangenome"
    Int cpu = 8
    Int memory_gb = 16
    Int heatmap_width = 1200
    Int heatmap_height = 800
  }

  command <<<
    set -euo pipefail

    # Always provision output dir and the files Cromwell expects
    mkdir -p final_output
    : > final_output/gene_presence_absence.csv
    : > final_output/summary_statistics.txt
    : > final_output/accessory_binary_genes.fa
    : > final_output/core_gene_alignment.aln
    : > final_output/number_of_new_genes.Rtab
    : > final_output/number_of_conserved_genes.Rtab

    # Minimal HTML created up-front so it's *always* there
    cat > final_output/pangenome_report.html <<'HTML'
<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>Pangenome Analysis Report</title>
<style>body{font-family:Arial,Helvetica,sans-serif;margin:20px}h1{color:#2c3e50}</style>
</head>
<body>
  <h1>Pangenome Analysis</h1>
  <p id="status">Initialized. (This page will show Roary results if available.)</p>
  <img id="heatmap" src="gene_presence_heatmap.png" alt="Gene presence/absence heatmap" onerror="this.style.display='none'">
</body>
</html>
HTML

    # Create placeholder heatmap (with Pillow if available)
    if python3 -c "from PIL import Image" 2>/dev/null; then
      python3 <<EOF
from PIL import Image, ImageDraw
img = Image.new('RGB', (~{heatmap_width}, ~{heatmap_height}), color='white')
d = ImageDraw.Draw(img)
d.rectangle([(50, 50), (~{heatmap_width}-50, ~{heatmap_height}-50)], outline='gray')
d.text((~{heatmap_width}//2, ~{heatmap_height}//2), "No Data Available", fill='black', anchor='mm', align='center')
img.save('final_output/gene_presence_heatmap.png')
EOF
    else
      echo "Pillow not available, using fallback heatmap generation." >&2
      convert -size ~{heatmap_width}x~{heatmap_height} xc:white -fill gray -draw "rectangle 50,50 $((~{heatmap_width}-50)),$((~{heatmap_height}-50))" -fill black -pointsize 24 -gravity center -annotate 0 "No Data Available" final_output/gene_presence_heatmap.png 2>/dev/null || : > final_output/gene_presence_heatmap.png
    fi

    # If user disabled pangenome, finish with the minimal report
    if [ "~{do_pangenome}" != "true" ]; then
      sed -i 's/Initialized./Pangenome step was disabled by parameters./' final_output/pangenome_report.html || true
      exit 0
    fi

    # Turn the WDL array into a bash array
    gffs=( ~{sep=' ' annotation_input} )
    # If no GFFs, keep the minimal report and exit 0
    if [ "${#gffs[@]}" -eq 0 ] || [ -z "${gffs[0]// }" ]; then
      sed -i 's/Initialized./No annotation GFFs were provided, so pangenome was not run./' final_output/pangenome_report.html || true
      exit 0
    fi

    # Copy inputs to a local folder (validate they exist & non-empty)
    mkdir -p gff_inputs
    valid=0
    for g in "${gffs[@]}"; do
      if [ -s "$g" ]; then
        cp "$g" gff_inputs/
        valid=$((valid+1))
      fi
    done
    if [ "$valid" -eq 0 ]; then
      sed -i 's/Initialized./No valid (non-empty) GFFs found, skipping pangenome./' final_output/pangenome_report.html || true
      exit 0
    fi

    # Try Roary; soft-fail so the pre-created outputs remain
    if command -v roary >/dev/null 2>&1; then
      set +e
      roary -f ~{output_prefix}_results \
            -p ~{cpu} -e -n -v \
            -i 90 -cd 99 \
            gff_inputs/*.gff
      roary_rc=$?
      set -e
      if [ "$roary_rc" -eq 0 ] && [ -d "~{output_prefix}_results" ]; then
        outdir="~{output_prefix}_results"
        cp -f "$outdir"/gene_presence_absence.csv                final_output/ 2>/dev/null || true
        cp -f "$outdir"/summary_statistics.txt                  final_output/ 2>/dev/null || true
        cp -f "$outdir"/accessory_binary_genes.fa               final_output/ 2>/dev/null || true
        # core alignment can be nested; find it
        ca=$(find "$outdir" -name "core_gene_alignment.aln" | head -1 || true)
        [ -n "${ca:-}" ] && cp -f "$ca" final_output/core_gene_alignment.aln || true
        cp -f "$outdir"/number_of_new_genes.Rtab                final_output/ 2>/dev/null || true
        cp -f "$outdir"/number_of_conserved_genes.Rtab          final_output/ 2>/dev/null || true

        # Improved heatmap generation with better error handling and memory optimization
        if command -v Rscript >/dev/null 2>&1; then
          cat > final_output/generate_visualizations.R <<'RS'
# Enhanced heatmap generation with validation and memory optimization
tryCatch({
  # Set memory limits and force garbage collection
  options(future.globals.maxSize = ~{memory_gb} * 1024^3)
  gc()

  # Load required packages
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap package not available")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("RColorBrewer package not available")
  }

  # Read and validate data
  if (!file.exists("gene_presence_absence.csv")) {
    stop("Input file not found")
  }

  # Read data with reduced memory footprint
  df <- read.csv("gene_presence_absence.csv", header=TRUE,
                check.names=FALSE, stringsAsFactors=FALSE,
                colClasses = "character")

  # Identify sample columns (heuristic)
  non_sample_cols <- c("Gene", "Annotation", "No..isolates", "No..sequences",
                      "Avg.sequences.per.isolate", "Genome.Fragment",
                      "Order.within.Fragment", "Accessory.Fragment",
                      "Accessory.Order.with.Fragment", "QC", "Min..sequence.length",
                      "Max..sequence.length", "Avg..sequence.length")

  sample_cols <- setdiff(colnames(df), non_sample_cols)

  # If heuristic failed, use all columns except first
  if (length(sample_cols) < 2) {
    sample_cols <- colnames(df)[-1]
  }

  # Convert to binary matrix with subsampling for large datasets
  mat <- as.matrix(df[, sample_cols, drop=FALSE])
  mat[mat != ""] <- 1
  mat[mat == ""] <- 0
  mode(mat) <- "numeric"

  # Subsample if matrix is too large
  if (nrow(mat) > 10000) {
    set.seed(123)
    mat <- mat[sample(1:nrow(mat), 10000), ]
  }

  # Only plot if we have at least 2 samples and 2 genes
  if (ncol(mat) >= 2 && nrow(mat) >= 2) {
    # Use sparse matrices if available
    if (requireNamespace("Matrix", quietly = TRUE)) {
      mat <- Matrix::Matrix(mat, sparse = TRUE)
    }

    png("gene_presence_heatmap.png", width=~{heatmap_width}, height=~{heatmap_height})
    pheatmap::pheatmap(
      mat,
      color = RColorBrewer::brewer.pal(9, "Blues")[c(1,9)],
      cluster_rows=TRUE,
      cluster_cols=TRUE,
      show_rownames=FALSE,
      main="Gene Presence/Absence",
      fontsize_col = max(6, min(12, 300/ncol(mat))) # Dynamic font sizing
    )
    dev.off()
  } else {
    stop("Insufficient data for heatmap (need ≥2 samples and ≥2 genes)")
  }
}, error = function(e) {
  # Create error image if heatmap fails
  png("gene_presence_heatmap.png", width=~{heatmap_width}, height=~{heatmap_height})
  plot.new()
  text(0.5, 0.5, paste("Heatmap Error:", e$message), adj=c(0.5,0.5))
  dev.off()
})
RS
          ( cd final_output && Rscript generate_visualizations.R ) || true
        fi

        # Stamp success into the HTML
        sed -i 's/Initialized./Roary completed successfully./' final_output/pangenome_report.html || true
      else
        sed -i 's/Initialized./Roary was invoked but did not complete successfully; placeholders shown./' final_output/pangenome_report.html || true
      fi  # FIXED: Added this missing fi
    else
      sed -i 's/Initialized./Roary not available in container; placeholders shown./' final_output/pangenome_report.html || true
    fi
  >>>

  runtime {
    docker: "gmboowa/roary-pillow:0.4"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    maxRetries: 2
    preemptible: 2
    continueOnReturnCode: true
    timeout: "24 hours"
  }

  output {
    File gene_presence_absence = "final_output/gene_presence_absence.csv"
    File summary_statistics    = "final_output/summary_statistics.txt"
    File pangenome_report      = "final_output/pangenome_report.html"
    File accessory_binary      = "final_output/accessory_binary_genes.fa"
    File core_alignment        = "final_output/core_gene_alignment.aln"
    File gene_heatmap          = "final_output/gene_presence_heatmap.png"
    File new_genes_data        = "final_output/number_of_new_genes.Rtab"
    File conserved_genes_data  = "final_output/number_of_conserved_genes.Rtab"
  }
}
task CORE_PHYLOGENY {
  input {
    File? alignment
    Boolean do_phylogeny = true
    String tree_prefix = "phylogeny"
    String model = "-nt -gtr"
    Int cpu = 16
    Int memory_gb = 16
    Int bootstrap_replicates = 100
    Int disk_gb = 100
  }

  command <<<
    set -euo pipefail

    # Initialize output directory and log file
    mkdir -p phylogeny_results
    echo "Core Phylogeny Analysis Log - $(date)" > phylogeny.log
    echo "====================================" >> phylogeny.log
    echo "Runtime Parameters:" >> phylogeny.log
    echo "- do_phylogeny: ~{do_phylogeny}" >> phylogeny.log
    echo "- tree_prefix: ~{tree_prefix}" >> phylogeny.log
    echo "- model: ~{model}" >> phylogeny.log
    echo "- cpu: ~{cpu}" >> phylogeny.log
    echo "- memory_gb: ~{memory_gb}" >> phylogeny.log
    echo "- bootstrap_replicates: ~{bootstrap_replicates}" >> phylogeny.log
    echo "====================================" >> phylogeny.log

    # Skip condition 1: User explicitly disabled phylogeny
    if [ "~{do_phylogeny}" != "true" ]; then
      echo "Phylogeny analysis disabled by user parameter" >> phylogeny.log
      echo "(PHYLOGENY_DISABLED);" > "phylogeny_results/core_~{tree_prefix}.nwk"
      echo "Analysis skipped by user request" > "phylogeny_results/core_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      exit 0
    fi

    # Skip condition 2: Alignment file not provided
    if [ ! -f "~{alignment}" ]; then
      echo "ERROR: Alignment file not found at path: ~{alignment}" >> phylogeny.log
      echo "(MISSING_ALIGNMENT);" > "phylogeny_results/core_~{tree_prefix}.nwk"
      echo "Alignment file missing" > "phylogeny_results/core_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      exit 0
    fi

    # Skip condition 3: Alignment file exists but is empty
    if [ ! -s "~{alignment}" ]; then
      echo "ERROR: Alignment file is empty: ~{alignment}" >> phylogeny.log
      echo "(EMPTY_ALIGNMENT);" > "phylogeny_results/core_~{tree_prefix}.nwk"
      echo "Alignment file empty" > "phylogeny_results/core_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      exit 0
    fi

    # Validate alignment content
    echo "Validating alignment file..." >> phylogeny.log
    seq_count=$(grep -c '^>' "~{alignment}" || echo 0)
    echo "Found $seq_count sequences in alignment" >> phylogeny.log

    # Skip condition 4: Insufficient sequences
    if [ "$seq_count" -lt 4 ]; then
      echo "ERROR: Insufficient sequences ($seq_count) for phylogenetic analysis (minimum 4 required)" >> phylogeny.log
      echo "(INSUFFICIENT_SEQUENCES_$seq_count);" > "phylogeny_results/core_~{tree_prefix}.nwk"
      echo "Only $seq_count sequences found" > "phylogeny_results/core_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      exit 0
    fi

    # Run phylogenetic analysis
    echo "Starting phylogenetic analysis with FastTree..." >> phylogeny.log
    echo "System memory information:" >> phylogeny.log
    free -h >> phylogeny.log

    set +e
    ulimit -v $((~{memory_gb} * 1024 * 1024))

    echo "Command: FastTree ~{model} -gamma -quiet -boot ~{bootstrap_replicates} \\" >> phylogeny.log
    echo "  -log phylogeny_results/core_~{tree_prefix}.log \\" >> phylogeny.log
    echo "  < ~{alignment} > phylogeny_results/core_~{tree_prefix}.nwk" >> phylogeny.log

    FastTree ~{model} \
      -gamma \
      -quiet \
      -boot ~{bootstrap_replicates} \
      -log "phylogeny_results/core_~{tree_prefix}.log" \
      < "~{alignment}" > "phylogeny_results/core_~{tree_prefix}.nwk" 2>> phylogeny.log
    exit_code=$?
    set -e

    # Handle FastTree failure
    if [ $exit_code -ne 0 ]; then
      echo "WARNING: FastTree exited with code $exit_code" >> phylogeny.log

      # Attempt reduced bootstrap replicates if memory issue
      if grep -qi "oom" phylogeny.log || grep -qi "killed" phylogeny.log; then
        echo "Attempting with reduced bootstrap replicates (50)..." >> phylogeny.log
        FastTree ~{model} \
          -gamma \
          -quiet \
          -boot 50 \
          -log "phylogeny_results/core_~{tree_prefix}_reduced.log" \
          < "~{alignment}" > "phylogeny_results/core_~{tree_prefix}.nwk" 2>> phylogeny.log || {
            echo "Fallback failed, generating minimal tree" >> phylogeny.log
            echo "(FASTTREE_FAILED);" > "phylogeny_results/core_~{tree_prefix}.nwk"
          }
      else
        echo "Generating minimal tree after failure" >> phylogeny.log
        echo "(FASTTREE_FAILED);" > "phylogeny_results/core_~{tree_prefix}.nwk"
      fi
    fi

    # Final validation of output
    if [ ! -s "phylogeny_results/core_~{tree_prefix}.nwk" ]; then
      echo "ERROR: Tree file is empty, generating minimal tree" >> phylogeny.log
      echo "(EMPTY_OUTPUT);" > "phylogeny_results/core_~{tree_prefix}.nwk"
    fi

    echo "Phylogenetic analysis completed at $(date)" >> phylogeny.log
    echo "Final output files:" >> phylogeny.log
    ls -lh phylogeny_results/* >> phylogeny.log
  >>>

  runtime {
    docker: "staphb/fasttree:2.1.11"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
    continueOnReturnCode: true
    timeout: "24 hours"
  }

  output {
    File phylogeny_tree = if (defined("phylogeny_results/core_~{tree_prefix}.nwk") &&
                           size("phylogeny_results/core_~{tree_prefix}.nwk") > 0)
                         then "phylogeny_results/core_~{tree_prefix}.nwk"
                         else write_lines(["();"])
    File? phylogeny_log_reduced = "phylogeny_results/core_~{tree_prefix}_reduced.log"
    File? phylogeny_log = "phylogeny_results/core_~{tree_prefix}.log"
    File execution_log = "phylogeny.log"
  }
}
task ACCESSORY_PHYLOGENY {
  input {
    File? alignment
    Boolean do_phylogeny = true
    String tree_prefix = "phylogeny"
    String model = "-nt -gtr"
    Int cpu = 8
    Int bootstrap_replicates = 100
    Int memory_gb = 16
  }

  command <<<
    #!/bin/bash
    set -euo pipefail

    # Initialize directories with proper permissions
    mkdir -p phylogeny_results
    chmod 777 phylogeny_results  # Ensure Docker can write

    # Start logging
    {
      echo "Accessory Phylogeny Analysis Log - $(date)"
      echo "========================================="
      echo "Runtime Parameters:"
      echo "- do_phylogeny: ~{do_phylogeny}"
      echo "- tree_prefix: ~{tree_prefix}"
      echo "- model: ~{model}"
      echo "- cpu: ~{cpu}"
      echo "- bootstrap_replicates: ~{bootstrap_replicates}"
      echo "- memory_gb: ~{memory_gb}"
      echo "========================================="
      echo "System Info:"
      free -h
      echo "-----------------------------------------"
    } > phylogeny.log

    # Skip conditions
    if [ "~{do_phylogeny}" = "false" ]; then
      echo "Accessory phylogeny disabled by user parameter" >> phylogeny.log
      echo "(ACCESSORY_PHYLOGENY_DISABLED);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Analysis skipped by user request" > "phylogeny_results/accessory_~{tree_prefix}.log"
      exit 0
    fi

    if [ -z "~{alignment}" ] || [ ! -f "~{alignment}" ]; then
      echo "ERROR: Accessory alignment file not found" >> phylogeny.log
      echo "(MISSING_ACCESSORY_ALIGNMENT);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Accessory alignment file missing" > "phylogeny_results/accessory_~{tree_prefix}.log"
      exit 0
    fi

    if [ ! -s "~{alignment}" ]; then
      echo "ERROR: Accessory alignment file is empty" >> phylogeny.log
      echo "(EMPTY_ACCESSORY_ALIGNMENT);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Accessory alignment file empty" > "phylogeny_results/accessory_~{tree_prefix}.log"
      exit 0
    fi

    # Validate alignment
    seq_count=$(grep -c '^>' "~{alignment}" || echo 0)
    echo "Found $seq_count sequences in accessory alignment" >> phylogeny.log

    if [ "$seq_count" -lt 4 ]; then
      echo "ERROR: Insufficient sequences ($seq_count)" >> phylogeny.log
      echo "(INSUFFICIENT_ACCESSORY_SEQUENCES_$seq_count);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Only $seq_count accessory sequences found" > "phylogeny_results/accessory_~{tree_prefix}.log"
      exit 0
    fi

    # Main analysis
    echo "Starting FastTree with ~{bootstrap_replicates} replicates..." >> phylogeny.log
    ulimit -v $((~{memory_gb} * 1024 * 1024))

    FastTree ~{model} \
      -gamma \
      -quiet \
      -boot ~{bootstrap_replicates} \
      -log "phylogeny_results/accessory_~{tree_prefix}.log" \
      < "~{alignment}" > "phylogeny_results/accessory_~{tree_prefix}.nwk" 2> "phylogeny_results/accessory_error.log" || {

      # Fallback for memory issues
      echo "FastTree failed, attempting with reduced replicates..." >> phylogeny.log
      FastTree ~{model} \
        -gamma \
        -quiet \
        -boot 50 \
        < "~{alignment}" > "phylogeny_results/accessory_~{tree_prefix}.nwk" 2>> "phylogeny.log" || {
          echo "Fallback failed, generating minimal tree" >> phylogeny.log
          echo "(ACCESSORY_TREE_FAILED);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
        }
    }

    # Final validation
    if [ ! -s "phylogeny_results/accessory_~{tree_prefix}.nwk" ]; then
      echo "ERROR: Output tree is empty" >> phylogeny.log
      echo "(EMPTY_ACCESSORY_OUTPUT);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
    fi

    echo "Analysis completed at $(date)" >> phylogeny.log
  >>>

  runtime {
    docker: "staphb/fasttree:2.1.11"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    preemptible: 2  # Allow for preemption
    maxRetries: 2   # Retry once on failure
  }

  output {
    File phylogeny_tree = "phylogeny_results/accessory_~{tree_prefix}.nwk"
    File phylogeny_log = if defined(glob("phylogeny_results/accessory_~{tree_prefix}.log"))
                         then "phylogeny_results/accessory_~{tree_prefix}.log"
                         else "phylogeny.log"
    File error_log = "phylogeny_results/accessory_error.log"
    File execution_log = "phylogeny.log"
  }
}

task MLST {
  input {
    Array[File]? assembly_output
    Boolean do_mlst = true
    Int cpu = 2
  }

  command <<<
    set -euo pipefail
    set -x

    # Initialize directories
    mkdir -p input_files mlst_results html_results

    # Start logging
    echo "MLST Analysis Log - $(date)" > mlst.log
    echo "==========================" >> mlst.log
    echo "Runtime Parameters:" >> mlst.log
    echo "- do_mlst: ~{do_mlst}" >> mlst.log
    echo "- cpu: ~{cpu}" >> mlst.log
    echo "==========================" >> mlst.log

    # Skip condition 1: User explicitly disabled MLST
    if [ "~{do_mlst}" != "true" ]; then
      echo "MLST analysis disabled by user parameter" >> mlst.log
      echo "MLST skipped by user request" > mlst_results/skipped.txt
      echo "<h1>MLST analysis skipped by user request</h1>" > html_results/skipped.html
      exit 0
    fi

    # Skip condition 2: No input files provided
    if [ -z "~{sep=' ' assembly_output}" ]; then
      echo "WARNING: No assembly files provided for MLST analysis" >> mlst.log
      echo "NO_INPUT_FILES" > mlst_results/skipped.txt
      echo "<h1>MLST analysis skipped - no input files provided</h1>" > html_results/skipped.html
      exit 0
    fi

    # Process input files
    echo "Input files verification:" >> mlst.log
    missing_files=0
    valid_files=0

    for f in ~{sep=' ' assembly_output}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> mlst.log
        missing_files=1
      else
        size=$(wc -c < "$f")
        echo "- Found input file: $f ($size bytes)" >> mlst.log
        if [ $size -gt 0 ]; then
          cp "$f" input_files/
          valid_files=$((valid_files + 1))
        else
          echo "WARNING: Empty input file: $f" >> mlst.log
        fi
      fi
    done

    # Skip condition 3: All input files missing or empty
    if [ $valid_files -eq 0 ]; then
      echo "ERROR: No valid input files provided" >> mlst.log
      echo "NO_VALID_INPUT_FILES" > mlst_results/skipped.txt
      echo "<h1>MLST analysis skipped - no valid input files</h1>" > html_results/skipped.html
      exit 0
    fi

    # Skip condition 4: Some files missing but others valid
    if [ $missing_files -ne 0 ]; then
      echo "WARNING: Some input files were missing, proceeding with available files" >> mlst.log
    fi

    # HTML template components
    HTML_HEADER='<!DOCTYPE html>
<html>
<head>
    <title>MLST Results</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .st { font-weight: bold; color: #e74c3c; }
        .warning { background-color: #fff3cd; }
        .error { background-color: #f8d7da; }
    </style>
</head>
<body>
    <h1>MLST Results</h1>
    <table>'

    HTML_FOOTER="    </table>
    <p>Generated at $(date)</p>
</body>
</html>"

    # Process each assembly file
    processed_files=0
    for asm_file in input_files/*; do
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//')
      output_file="mlst_results/${sample_name}_mlst.tsv"
      debug_file="mlst_results/${sample_name}_debug.log"
      temp_file="mlst_results/${sample_name}_temp.tsv"

      echo "Processing $sample_name" >> mlst.log

      if mlst --debug "$asm_file" > "$temp_file" 2> "$debug_file"; then
        {
          echo -e "sampleID\tspecies\tST\tgapA\tinfB\tmdh\tpgi\tphoE\trpoB\ttonB"
          awk -v sample="$sample_name" '
function extract_allele(field,   a) {
  split(field, a, /[()]/);
  return (a[2] != "") ? a[2] : "N/A";
}
BEGIN { OFS = "\t" }
NR==1 && $1=="FILE" { next }
NF > 0 {
  species = $2;
  st = $3; sub(/^ST/, "", st);
  gapA = extract_allele($4);
  infB = extract_allele($5);
  mdh  = extract_allele($6);
  pgi  = extract_allele($7);
  phoE = extract_allele($8);
  rpoB = extract_allele($9);
  tonB = extract_allele($10);
  print sample, species, st, gapA, infB, mdh, pgi, phoE, rpoB, tonB;
}' "$temp_file"
        } > "$output_file"

        if [ ! -s "$output_file" ] || [ "$(wc -l < "$output_file")" -le 1 ]; then
          echo "WARNING: MLST output processing failed for $asm_file" >> mlst.log
          echo -e "sampleID\tspecies\tST\tgapA\tinfB\tmdh\tpgi\tphoE\trpoB\ttonB" > "$output_file"
          echo -e "$sample_name\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR" >> "$output_file"
        else
          processed_files=$((processed_files + 1))
        fi
      else
        echo "WARNING: MLST failed for $asm_file" >> mlst.log
        echo "Debug output:" >> mlst.log
        cat "$debug_file" >> mlst.log
        echo -e "sampleID\tspecies\tST\tgapA\tinfB\tmdh\tpgi\tphoE\trpoB\ttonB" > "$output_file"
        echo -e "$sample_name\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR" >> "$output_file"
      fi

      [ -f "$temp_file" ] && rm -f "$temp_file"
    done

    # Generate HTML reports
    for tsv_file in mlst_results/*_mlst.tsv; do
      sample_name=$(basename "$tsv_file" | sed 's/_mlst.tsv//')
      html_file="html_results/${sample_name}_mlst.html"

      {
        echo "$HTML_HEADER"
        echo "<thead><tr>"
        echo "<th>Sample ID</th>"
        echo "<th>Species</th>"
        echo "<th>ST</th>"
        echo "<th>gapA</th>"
        echo "<th>infB</th>"
        echo "<th>mdh</th>"
        echo "<th>pgi</th>"
        echo "<th>phoE</th>"
        echo "<th>rpoB</th>"
        echo "<th>tonB</th>"
        echo "</tr></thead><tbody>"

        line_count=$(wc -l < "$tsv_file")
        if [ $line_count -gt 1 ]; then
          tail -n +2 "$tsv_file" | while IFS=$'\t' read -r sample species st gapA infB mdh pgi phoE rpoB tonB; do
            echo "<tr>"
            echo "<td>$sample</td>"
            echo "<td>$species</td>"
            if [[ "$st" == "ERROR" ]]; then
              echo "<td class=\"error\">$st</td>"
            else
              echo "<td class=\"st\">$st</td>"
            fi
            echo "<td>$gapA</td>"
            echo "<td>$infB</td>"
            echo "<td>$mdh</td>"
            echo "<td>$pgi</td>"
            echo "<td>$phoE</td>"
            echo "<td>$rpoB</td>"
            echo "<td>$tonB</td>"
            echo "</tr>"
          done
        else
          echo "<tr class='warning'>"
          echo "<td colspan='10'>No MLST results for $sample_name</td>"
          echo "</tr>"
        fi

        echo "</tbody>"
        echo "$HTML_FOOTER"
      } > "$html_file"
    done

    # Create combined results
    echo -e "sampleID\tspecies\tST\tgapA\tinfB\tmdh\tpgi\tphoE\trpoB\ttonB" > mlst_results/combined_mlst.tsv
    for tsv_file in mlst_results/*_mlst.tsv; do
      if [ -s "$tsv_file" ] && [ $(wc -l < "$tsv_file") -gt 1 ]; then
        tail -n +2 "$tsv_file" >> mlst_results/combined_mlst.tsv
      fi
    done

    # De-duplicate combined results
    awk 'NR==1{print; next} !seen[$0]++' mlst_results/combined_mlst.tsv > mlst_results/combined_mlst.dedup.tsv && \
    mv mlst_results/combined_mlst.dedup.tsv mlst_results/combined_mlst.tsv

    # Generate combined HTML
    {
      echo "$HTML_HEADER"
      echo "<thead><tr>"
      echo "<th>Sample ID</th>"
      echo "<th>Species</th>"
      echo "<th>ST</th>"
      echo "<th>gapA</th>"
      echo "<th>infB</th>"
      echo "<th>mdh</th>"
      echo "<th>pgi</th>"
      echo "<th>phoE</th>"
      echo "<th>rpoB</th>"
      echo "<th>tonB</th>"
      echo "</tr></thead><tbody>"

      if [ -s "mlst_results/combined_mlst.tsv" ]; then
        while IFS=$'\t' read -r sample species st gapA infB mdh pgi phoE rpoB tonB; do
          if [ "$sample" == "sampleID" ]; then
            continue
          fi
          echo "<tr>"
          echo "<td>$sample</td>"
          echo "<td>$species</td>"
          if [[ "$st" == "ERROR" ]]; then
            echo "<td class=\"error\">$st</td>"
          else
            echo "<td class=\"st\">$st</td>"
          fi
          echo "<td>$gapA</td>"
          echo "<td>$infB</td>"
          echo "<td>$mdh</td>"
          echo "<td>$pgi</td>"
          echo "<td>$phoE</td>"
          echo "<td>$rpoB</td>"
          echo "<td>$tonB</td>"
          echo "</tr>"
        done < mlst_results/combined_mlst.tsv
      else
        echo "<tr class='error'>"
        echo "<td colspan='10'>No valid MLST results were generated</td>"
        echo "</tr>"
      fi

      echo "</tbody>"
      echo "$HTML_FOOTER"
    } > html_results/combined_mlst.html

    echo "MLST completed at $(date)" >> mlst.log
    echo "Processed $processed_files files successfully" >> mlst.log
    echo "Output files:" >> mlst.log
    ls -lh mlst_results/* html_results/* >> mlst.log 2>&1 || true
  >>>

  runtime {
    docker: "staphb/mlst:2.19.0"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
    timeout: "6 hours"
  }

  output {
    Array[File] mlst_outputs = if do_mlst then glob("mlst_results/*_mlst.tsv") else ["mlst_results/skipped.txt"]
    File combined_mlst = if do_mlst then "mlst_results/combined_mlst.tsv" else "mlst_results/skipped.txt"
    File combined_html = if do_mlst then "html_results/combined_mlst.html" else "html_results/skipped.html"
    File mlst_log = "mlst.log"
    Array[File] debug_logs = if do_mlst then glob("mlst_results/*_debug.log") else []
  }
}
task VARIANT_CALLING {
  input {
    Array[File]+ input_reads
    File reference_genome
    Boolean do_variant_calling = true
    String reference_type = "genbank"
    Int cpu = 8
    Int min_quality = 20
    Int memory_gb = 8
    Int max_retries = 2
  }

  command <<<
    #!/bin/bash
    set -euo pipefail
    set -x

    # Initialize directories and logs
    mkdir -p variants || { echo "ERROR: Failed to create variants directory" >&2; exit 1; }
    echo "Variant Calling Analysis Log - $(date)" > variant.log
    echo "=====================================" >> variant.log
    echo "Runtime Parameters:" >> variant.log
    echo "- do_variant_calling: ~{do_variant_calling}" >> variant.log
    echo "- reference_type: ~{reference_type}" >> variant.log
    echo "- cpu: ~{cpu}" >> variant.log
    echo "- min_quality: ~{min_quality}" >> variant.log
    echo "- memory_gb: ~{memory_gb}" >> variant.log
    echo "=====================================" >> variant.log

    # Function to ensure HTML output exists
    ensure_output() {
      local status=$1
      mkdir -p variants

      if [ "$status" == "success" ]; then
        if [ ! -f variants/variant_summary.html ]; then
          cat > variants/variant_summary.html <<EOF
<!DOCTYPE html>
<html>
<head><title>Variant Calling Summary</title></head>
<body>
  <h1>Variant Calling Summary</h1>
  <p>Analysis completed but final report not generated</p>
  <p>See variant.log for details</p>
</body>
</html>
EOF
        fi
      else
        cat > variants/variant_summary.html <<EOF
<!DOCTYPE html>
<html>
<head><title>Variant Calling Summary</title></head>
<body>
  <h1>Variant Calling Summary</h1>
  <p>Analysis failed - see variant.log for details</p>
</body>
</html>
EOF
      fi
      chmod 644 variants/variant_summary.html
    }

    # Handle skip conditions
    if [ "~{do_variant_calling}" != "true" ]; then
      echo "Variant calling disabled by user parameter" >> variant.log
      ensure_output "skipped"
      echo "Variant calling skipped by user request" > variants/skipped.txt
      exit 0
    fi

    # Validate inputs
    echo "Validating input files..." >> variant.log
    valid_files=0
    files=()
    for f in ~{sep=' ' input_reads}; do
      if [ ! -f "$f" ]; then
        echo "ERROR: Input file not found: $f" >> variant.log
        exit 1
      fi
      size=$(wc -c < "$f")
      if [ $size -gt 0 ]; then
        files+=("$f")
        valid_files=$((valid_files + 1))
        echo "Valid input: $f ($size bytes)" >> variant.log
      else
        echo "ERROR: Empty input file: $f" >> variant.log
        exit 1
      fi
    done

    if [ $valid_files -eq 0 ]; then
      echo "ERROR: No valid input files" >> variant.log
      ensure_output "failed"
      exit 1
    fi

    if [ $((valid_files % 2)) -ne 0 ]; then
      echo "ERROR: Odd number of input files (needs pairs)" >> variant.log
      ensure_output "failed"
      exit 1
    fi

    if [ ! -f "~{reference_genome}" ]; then
      echo "ERROR: Reference genome not found" >> variant.log
      ensure_output "failed"
      exit 1
    fi

    # Create initial empty HTML as placeholder
    ensure_output "running"

    # Process samples with retry logic
    RESULTS_TEMP=$(mktemp) || { echo "ERROR: Failed to create temp file" >> variant.log; exit 1; }
    echo -e "Sample\tStatus\tSNPs\tINDELs\tTotal\tAttempts" > "$RESULTS_TEMP"
    processed_samples=0

    for ((i=0; i<${#files[@]}; i+=2)); do
      R1="${files[i]}"
      R2="${files[i+1]}"
      SAMPLE_NAME=$(basename "$R1" | sed 's/[._][Rr]1.*//; s/[._]1.*//; s/[._][12].*//')
      SAMPLE_DIR="variants/$SAMPLE_NAME"

      mkdir -p "$SAMPLE_DIR" || {
        echo "ERROR: Failed to create directory for $SAMPLE_NAME" >> variant.log
        echo -e "$SAMPLE_NAME\tfailed\t0\t0\t0\t0" >> "$RESULTS_TEMP"
        continue
      }

      echo "Processing $SAMPLE_NAME" >> variant.log
      SNP_COUNT=0
      INDEL_COUNT=0
      TOTAL=0
      STATUS="failed"  # Default to failed unless proven otherwise
      ATTEMPTS=0
      MAX_ATTEMPTS=~{max_retries}

      # Handle reference type
      if [ "~{reference_type}" == "genbank" ]; then
        ref_path="$SAMPLE_DIR/reference.gbk"
      else
        ref_path="$SAMPLE_DIR/reference.fasta"
      fi

      if ! cp "~{reference_genome}" "$ref_path"; then
        echo "ERROR: Failed to copy reference for $SAMPLE_NAME" >> variant.log
        echo -e "$SAMPLE_NAME\tfailed\t0\t0\t0\t0" >> "$RESULTS_TEMP"
        continue
      fi

      # Run variant calling with retry logic
      while [ $ATTEMPTS -le $MAX_ATTEMPTS ] && [ "$STATUS" != "success" ]; do
        ATTEMPTS=$((ATTEMPTS + 1))
        echo "Attempt $ATTEMPTS for $SAMPLE_NAME" >> variant.log

        # Clean previous attempt files if they exist
        rm -f "$SAMPLE_DIR/$SAMPLE_NAME.bam" "$SAMPLE_DIR/$SAMPLE_NAME.bam.bai" "$SAMPLE_DIR/$SAMPLE_NAME.vcf" 2>/dev/null || true

        set +e
        snippy --cpus ~{cpu} --minqual ~{min_quality} \
               --ref "$ref_path" --R1 "$R1" --R2 "$R2" \
               --outdir "$SAMPLE_DIR" --prefix "$SAMPLE_NAME" --force 2>> variant.log
        snippy_exit=$?
        set -e

        if [ $snippy_exit -eq 0 ] && [ -s "$SAMPLE_DIR/$SAMPLE_NAME.vcf" ]; then
          STATUS="success"
          SNP_COUNT=$(grep -c 'TYPE=snp' "$SAMPLE_DIR/$SAMPLE_NAME.vcf" 2>/dev/null || echo 0)
          INDEL_COUNT=$(awk '/TYPE=indel/ {count++} END {print count+0}' "$SAMPLE_DIR/$SAMPLE_NAME.vcf" 2>/dev/null)
          TOTAL=$((SNP_COUNT + INDEL_COUNT))
          processed_samples=$((processed_samples + 1))
          echo "Successfully processed $SAMPLE_NAME on attempt $ATTEMPTS" >> variant.log
        else
          echo "WARNING: snippy failed for $SAMPLE_NAME on attempt $ATTEMPTS (exit $snippy_exit)" >> variant.log
          if [ $ATTEMPTS -lt $MAX_ATTEMPTS ]; then
            echo "Will retry $SAMPLE_NAME after a short delay..." >> variant.log
            sleep $((ATTEMPTS * 30))  # Exponential backoff
          fi
        fi
      done

      echo -e "$SAMPLE_NAME\t$STATUS\t$SNP_COUNT\t$INDEL_COUNT\t$TOTAL\t$ATTEMPTS" >> "$RESULTS_TEMP"
      echo "Processed $SAMPLE_NAME: Status=$STATUS, SNPs=$SNP_COUNT, INDELs=$INDEL_COUNT, Total=$TOTAL, Attempts=$ATTEMPTS" >> variant.log
    done

    # Generate final HTML report
    cat > variants/variant_summary.html <<EOF
<!DOCTYPE html>
<html>
<head>
  <title>Variant Calling Summary</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    table { border-collapse: collapse; width: 100%; margin-top: 20px; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }
    th { background-color: #3498db; color: white; }
    tr:nth-child(even) { background-color: #f2f2f2; }
    .success { color: #27ae60; }
    .warning { color: #f39c12; }
    .error { color: #e74c3c; }
  </style>
</head>
<body>
  <h1>Variant Calling Summary</h1>
  <p>Generated: $(date)</p>
  <p>Samples processed: $processed_samples</p>
  <table>
    <thead>
      <tr>
        <th>Sample</th>
        <th>Status</th>
        <th>SNPs</th>
        <th>INDELs</th>
        <th>Total Variants</th>
        <th>Attempts</th>
      </tr>
    </thead>
    <tbody>
EOF

    while IFS=$'\t' read -r sample status snps indels total attempts; do
      if [ "$sample" == "Sample" ]; then continue; fi
      if [ "$status" == "failed" ]; then
        echo "<tr class=\"error\"><td>$sample</td><td>Failed</td><td>-</td><td>-</td><td>-</td><td>$attempts</td></tr>" >> variants/variant_summary.html
      else
        if [ $total -eq 0 ]; then
          row_class="success"
        elif [ $total -lt 100 ]; then
          row_class="warning"
        else
          row_class="error"
        fi
        echo "<tr><td>$sample</td><td>Success</td><td>$snps</td><td>$indels</td><td class=\"$row_class\">$total</td><td>$attempts</td></tr>" >> variants/variant_summary.html
      fi
    done < "$RESULTS_TEMP"

    cat >> variants/variant_summary.html <<EOF
    </tbody>
  </table>
</body>
</html>
EOF

    # Final cleanup
    rm -f "$RESULTS_TEMP"
    echo "Variant calling completed at $(date)" >> variant.log
    echo "Output files:" >> variant.log
    find variants -type f -exec ls -lh {} \; >> variant.log 2>/dev/null || true

    # Exit with error if any samples failed
    if grep -q "failed" "$RESULTS_TEMP"; then
      echo "WARNING: Some samples failed processing" >> variant.log
      exit 1
    fi
  >>>

  runtime {
    docker: "quay.io/biocontainers/snippy:4.6.0--hdfd78af_1"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
    preemptible: 2
    timeout: "24 hours"
  }

  output {
    Array[File] vcf_files        = glob("variants/*/*.vcf")
    Array[File] variant_output   = glob("variants/*/*")
    Array[File] variant_dirs     = glob("variants/*")
    File        variant_log      = "variant.log"
    File        variant_summary  = "variants/variant_summary.html"
    String      variants_dir     = "variants"
  }
}
task AMR_PROFILING {
  input {
    Array[File]? assembly_output
    Boolean do_amr_profiling = true
    File? local_db
    Boolean use_local_db = false
    Int minid = 90
    Int mincov = 80
    Int cpu = 2
    Boolean abricate_nopath = true
    Boolean abricate_make_summary = false
    Int merge_minid = 95
    Int merge_mincov = 90  
  }

  command <<<
    set -euo pipefail
    set -x

    # Initialize directories
    mkdir -p amr_results html_results amr_results/raw_csv
    echo "AMR Profiling Analysis Log - $(date)" > amr.log
    echo "===================================" >> amr.log
    echo "Runtime Parameters:" >> amr.log
    echo "- do_amr_profiling: ~{do_amr_profiling}" >> amr.log
    echo "- use_local_db: ~{use_local_db}" >> amr.log
    echo "- minid: ~{minid}" >> amr.log
    echo "- mincov: ~{mincov}" >> amr.log
    echo "- cpu: ~{cpu}" >> amr.log
    echo "- abricate_nopath: ~{abricate_nopath}" >> amr.log
    echo "- abricate_make_summary: ~{abricate_make_summary}" >> amr.log
    echo "- merge_minid: ~{merge_minid}" >> amr.log
    echo "- merge_mincov: ~{merge_mincov}" >> amr.log
    echo "===================================" >> amr.log

    # HTML template
    HTML_HEADER='<!DOCTYPE html>
<html>
<head>
    <title>AMR Profiling Results</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; position: sticky; top: 0; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .gene { font-weight: bold; color: #e74c3c; }
        .antibiotic { font-style: italic; color: #27ae60; }
        .coverage-high { background-color: #e6ffe6; }
        .coverage-medium { background-color: #fff9e6; }
        .coverage-low { background-color: #ffe6e6; }
        .summary-card {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .error-row { background-color: #f8d7da; }
    </style>
</head>
<body>
    <h1>Antimicrobial Resistance Profiling Results</h1>
    <div class="summary-card">
        <p><strong>Analysis Date:</strong> '"$(date)"'</p>
        <p><strong>Database Used:</strong> ~{if use_local_db then "Custom" else "ResFinder"}</p>
        <p><strong>Parameters:</strong></p>
        <ul>
            <li>Minimum Identity: ~{minid}%</li>
            <li>Minimum Coverage: ~{mincov}%</li>
            <li>Threads: ~{cpu}</li>
        </ul>
    </div>
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Contig</th>
                <th>Gene</th>
                <th>%Coverage</th>
                <th>%Identity</th>
                <th>Product</th>
                <th>Resistance</th>
            </tr>
        </thead>
        <tbody>'

    HTML_FOOTER='        </tbody>
    </table>
</body>
</html>'

    # Skip conditions
    if [ "~{do_amr_profiling}" != "true" ]; then
      echo "AMR profiling disabled by user parameter" >> amr.log
      echo "AMR profiling skipped by user request" > amr_results/skipped.txt
      echo "<h1>AMR profiling skipped by user request</h1>" > html_results/skipped.html
      exit 0
    fi

    if [ -z "~{sep=' ' assembly_output}" ]; then
      echo "ERROR: No assembly files provided for AMR profiling" >> amr.log
      echo "NO_INPUT_FILES" > amr_results/skipped.txt
      echo "<h1>AMR profiling skipped - no input files provided</h1>" > html_results/skipped.html
      exit 0
    fi

    # Verify input files
    echo "Input files verification:" >> amr.log
    valid_files=0
    for f in ~{sep=' ' assembly_output}; do
      if [ ! -f "$f" ]; then
        echo "WARNING: Input file not found: $f" >> amr.log
      else
        size=$(wc -c < "$f")
        if [ $size -gt 0 ]; then
          echo "- Valid input file: $f ($size bytes)" >> amr.log
          valid_files=$((valid_files + 1))
        else
          echo "WARNING: Empty input file: $f" >> amr.log
        fi
      fi
    done

    if [ $valid_files -eq 0 ]; then
      echo "ERROR: No valid input files available" >> amr.log
      echo "NO_VALID_INPUTS" > amr_results/skipped.txt
      echo "<h1>AMR profiling skipped - no valid input files</h1>" > html_results/skipped.html
      exit 0
    fi

    # Database setup
    db_to_use="resfinder"
    if [ "~{use_local_db}" == "true" ]; then
      if [ ! -f "~{local_db}" ]; then
        echo "ERROR: Local database file not found" >> amr.log
        echo "MISSING_LOCAL_DB" > amr_results/skipped.txt
        echo "<h1>AMR profiling skipped - local database missing</h1>" > html_results/skipped.html
        exit 0
      fi
      mkdir -p /root/abricate/db/resfinder_db
      cp "~{local_db}" /root/abricate/db/resfinder_db/resfinder.fa || {
        echo "ERROR: Failed to copy local database" >> amr.log
        echo "DB_COPY_FAILED" > amr_results/skipped.txt
        echo "<h1>AMR profiling skipped - database setup failed</h1>" > html_results/skipped.html
        exit 0
      }
      abricate --setupdb --db resfinder --debug >> amr.log 2>&1 || {
        echo "ERROR: Failed to setup local AMR database" >> amr.log
        echo "DB_SETUP_FAILED" > amr_results/skipped.txt
        echo "<h1>AMR profiling skipped - database setup failed</h1>" > html_results/skipped.html
        exit 0
      }
    else
      abricate --list | grep -q "resfinder" || abricate --setupdb --db resfinder >> amr.log 2>&1 || {
        echo "ERROR: Failed to setup resfinder database" >> amr.log
        echo "DB_SETUP_FAILED" > amr_results/skipped.txt
        echo "<h1>AMR profiling skipped - database setup failed</h1>" > html_results/skipped.html
        exit 0
      }
    fi

    # Process samples (modified for consistent naming)
    processed_samples=0
    for asm_file in ~{sep=' ' assembly_output}; do
      [ ! -f "$asm_file" ] && continue
      [ ! -s "$asm_file" ] && continue

      sample_name=$(basename "$asm_file" | sed -E 's/\.(fa|fasta|fna|fsa|contigs|scaffolds)(\.(gz|bz2))?$//i')
      output_file="amr_results/${sample_name}_amr.tsv"
      html_file="html_results/${sample_name}_amr.html"
      raw_csv="amr_results/raw_csv/${sample_name}.csv"

      echo "Processing $sample_name" >> amr.log

      # Run abricate (unchanged)
      extra_flags=()
      if [ "~{abricate_nopath}" = "true" ]; then
        extra_flags+=("--nopath")
      fi
      abricate \
        --db $db_to_use \
        --minid ~{minid} \
        --mincov ~{mincov} \
        --threads ~{cpu} \
        --csv \
        "${extra_flags[@]}" \
        "$asm_file" > "${output_file}.tmp" 2>> amr.log || {
          echo "WARNING: AMR profiling failed for $sample_name" >> amr.log
          echo -e "SAMPLE\tCONTIG\tGENE\t%COVERAGE\t%IDENTITY\tPRODUCT\tRESISTANCE" > "$output_file"
          echo -e "$sample_name\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR" >> "$output_file"
          {
            echo "$HTML_HEADER"
            echo "<tr class=\"error-row\"><td colspan=\"7\">AMR profiling failed for $sample_name</td></tr>"
            echo "$HTML_FOOTER"
          } > "$html_file"
          continue
        }

      # Save raw CSV (unchanged)
      cp -f "${output_file}.tmp" "$raw_csv" || true

      # Convert CSV -> TSV (unchanged)
      awk -F',' -v sample="$sample_name" '
      NR==1 {
        for (i=1; i<=NF; i++) {
          h = $i; gsub(/^#*/,"",h); gsub(/"/,"",h); gsub(/^[[:space:]]+|[[:space:]]+$/, "", h); H[h]=i
        }
        OFS="\t"
        print "SAMPLE","CONTIG","GENE","%COVERAGE","%IDENTITY","PRODUCT","RESISTANCE"
        next
      }
      {
        contig = (H["SEQUENCE"] ? $(H["SEQUENCE"]) : "NA")
        gene   = (H["GENE"] ? $(H["GENE"]) : "NA")
        pcov   = (H["%COVERAGE"] ? $(H["%COVERAGE"]) : (H["COVERAGE"] ? $(H["COVERAGE"]) : "NA"))
        pid    = (H["%IDENTITY"] ? $(H["%IDENTITY"]) : "NA")
        prod   = (H["PRODUCT"] ? $(H["PRODUCT"]) : "NA")
        res    = (H["RESISTANCE"] ? $(H["RESISTANCE"]) : "N/A")
        gsub(/"/,"",contig); gsub(/"/,"",gene); gsub(/"/,"",pcov); gsub(/"/,"",pid); gsub(/"/,"",prod); gsub(/"/,"",res)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", contig)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", gene)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", pcov)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", pid)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", prod)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", res)
        print sample, contig, gene, pcov, pid, prod, res
      }' "${output_file}.tmp" > "$output_file"

      # --- FIX A: Normalize per-sample TSV to remove leading blanks / duplicate header ---
      sed -i '1{/^$/d;}' "$output_file"
      awk 'NR==1{print; h=$0; next} $0!=h{print}' "$output_file" > "${output_file}.clean" && mv "${output_file}.clean" "$output_file"
      # -------------------------------------------------------------------------------

      # Best-hit filtering (patched to guard against header-as-row)
      awk -F'\t' '
        NR==1 { header=$0; next }
        $1=="SAMPLE" && $3=="GENE" { next }  # FIX B: drop accidental header row
        {
          g=$3
          cov=$4; sub(/%$/,"",cov); if (cov=="") cov=0
          pid=$5; sub(/%$/,"",pid); if (pid=="") pid=0
          key=$1"\t"g
          if (!(key in best_cov) || cov>best_cov[key] || (cov==best_cov[key] && pid>best_id[key])) {
            best_cov[key]=cov
            best_id[key]=pid
            best_line[key]=$0
          }
        }
        END {
          print header
          for (k in best_line) print best_line[k]
        }
      ' "$output_file" | sort -t$'\t' -k1,1 -k3,3 > "${output_file}.best"
      mv -f "${output_file}.best" "$output_file"

      # Generate per-sample HTML
      {
        echo "$HTML_HEADER"
        tail -n +2 "$output_file" | while IFS=$'\t' read -r sample contig gene coverage identity product resistance; do
          # FIX C: skip if a stray header row sneaks through
          if [[ "$sample" == "SAMPLE" && "$gene" == "GENE" ]]; then
            continue
          fi
          if [[ "$gene" == *"ERROR"* ]]; then
            echo "<tr class=\"error-row\">"
            echo "<td>$sample</td>"
            echo "<td colspan=\"6\">AMR profiling failed</td>"
            echo "</tr>"
          else
            echo "<tr>"
            echo "<td>$sample</td>"
            echo "<td>$contig</td>"
            echo "<td class=\"gene\">$gene</td>"
            echo "<td>$coverage</td>"
            echo "<td>$identity</td>"
            echo "<td>$product</td>"
            echo "<td class=\"antibiotic\">$resistance</td>"
            echo "</tr>"
          fi
        done
        echo "$HTML_FOOTER"
      } > "$html_file"

      processed_samples=$((processed_samples + 1))
      rm -f "${output_file}.tmp" || true
    done

    # Optional abricate summary (unchanged)
    if [ "~{abricate_make_summary}" = "true" ]; then
      if ls amr_results/raw_csv/*.csv >/dev/null 2>&1; then
        abricate --summary amr_results/raw_csv/*.csv > amr_results/abricate_summary.tsv 2>> amr.log || true
        gzip -f amr_results/abricate_summary.tsv || true
      fi
    fi

    # Generate combined HTML report
    if [ $processed_samples -gt 0 ]; then
      {
        echo "$HTML_HEADER"
        # Read directly from per-sample TSVs (no combined_amr.tsv needed)
        for tsv_file in amr_results/*_amr.tsv; do
          tail -n +2 "$tsv_file" | while IFS=$'\t' read -r sample contig gene coverage identity product resistance; do
            # FIX C (also applied here): skip stray header rows
            if [[ "$sample" == "SAMPLE" && "$gene" == "GENE" ]]; then
              continue
            fi
            if [[ "$gene" == *"ERROR"* ]]; then
              echo "<tr class=\"error-row\">"
              echo "<td>$sample</td>"
              echo "<td colspan=\"6\">AMR profiling failed</td>"
              echo "</tr>"
            else
              # Apply merge thresholds here
              cov="${coverage//%}"
              id="${identity//%}"
              if (( cov >= ~{merge_mincov} )) && (( id >= ~{merge_minid} )); then
                echo "<tr>"
                echo "<td>$sample</td>"
                echo "<td>$contig</td>"
                echo "<td class=\"gene\">$gene</td>"
                echo "<td>$coverage</td>"
                echo "<td>$identity</td>"
                echo "<td>$product</td>"
                echo "<td class=\"antibiotic\">$resistance</td>"
                echo "</tr>"
              fi
            fi
          done
        done
        echo "$HTML_FOOTER"
      } > html_results/combined_amr.html
    else
      echo "ERROR: No samples processed successfully" >> amr.log
      echo "<h1>No samples processed successfully</h1>" > html_results/combined_amr.html
    fi

    echo "AMR profiling completed at $(date)" >> amr.log
    echo "Processed $processed_samples samples successfully" >> amr.log
  >>>

  runtime {
    docker: "staphb/abricate:1.0.0"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
    timeout: "10 hours"
  }

  output {
    # Per-sample outputs
    Array[File] html_reports = if do_amr_profiling then glob("html_results/*_amr.html")
                               else ["html_results/skipped.html"]
    Array[File] amr_outputs  = if do_amr_profiling then glob("amr_results/*_amr.tsv")
                               else ["amr_results/skipped.txt"]

    # Combined HTML
    File combined_html       = if do_amr_profiling then "html_results/combined_amr.html"
                               else "html_results/skipped.html"

    File amr_log             = "amr.log"

    # Backward-compatible alias so existing references to AMR_PROFILING.combined_amr keep working.
    File combined_amr        = combined_html
  }
}

task MGE {
  input {
    File assembly
    Boolean do_mge_analysis = true
    File? local_db
    Boolean use_local_db = false
    Int cpu = 4
  }

  command <<<
    set -euo pipefail
    set -x

    # Initialize directories
    mkdir -p mge_results html_results
    echo "Mobile Genetic Element Analysis Log - $(date)" > mge.log
    echo "===========================================" >> mge.log
    echo "Runtime Parameters:" >> mge.log
    echo "- do_mge_analysis: ~{do_mge_analysis}" >> mge.log
    echo "- use_local_db: ~{use_local_db}" >> mge.log
    echo "- cpu: ~{cpu}" >> mge.log
    echo "===========================================" >> mge.log

    # Skip condition: User explicitly disabled MGE analysis
    if [ "~{do_mge_analysis}" != "true" ]; then
      echo "MGE analysis disabled by user parameter" >> mge.log
      echo "MGE analysis skipped by user request" > mge_results/skipped.txt
      echo "<h1>MGE analysis skipped by user request</h1>" > html_results/skipped.html
      exit 0
    fi

    # Verify input file
    valid_files=0
    if [ ! -f "~{assembly}" ]; then
      echo "WARNING: Input file not found: ~{assembly}" >> mge.log
    else
      size=$(wc -c < "~{assembly}")
      if [ $size -gt 0 ]; then
        echo "- Valid input file: ~{assembly} ($size bytes)" >> mge.log
        valid_files=$((valid_files + 1))
      else
        echo "WARNING: Empty input file: ~{assembly}" >> mge.log
      fi
    fi

    # Skip condition: No valid input file
    if [ $valid_files -eq 0 ]; then
      echo "ERROR: No valid input file available" >> mge.log
      echo "NO_VALID_INPUTS" > mge_results/skipped.txt
      echo "<h1>MGE analysis skipped - no valid input file</h1>" > html_results/skipped.html
      exit 0
    fi

    # HTML template components (unchanged)
    HTML_HEADER='<!DOCTYPE html>
<html>
<head>
    <title>Mobile Genetic Element Analysis</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; position: sticky; top: 0; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .mge { font-weight: bold; color: #e74c3c; }
        .coverage-high { background-color: #e6ffe6; }
        .coverage-medium { background-color: #fff9e6; }
        .coverage-low { background-color: #ffe6e6; }
        .summary-card {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .error-row { background-color: #f8d7da; }
    </style>
</head>
<body>
    <h1>Mobile Genetic Element Analysis Report</h1>
    <div class="summary-card">
        <p><strong>Analysis Date:</strong> <span id="analysis-date"></span></p>
        <p><strong>Database Used:</strong> <span id="database-used"></span></p>
        <p><strong>Parameters:</strong></p>
        <ul>
            <li>Minimum Coverage: 80%</li>
            <li>Minimum Identity: 90%</li>
            <li>Threads: ~{cpu}</li>
        </ul>
        <p><strong>Samples Processed:</strong> <span id="samples-processed"></span></p>
    </div>
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>MGE Type</th>
                <th>Gene</th>
                <th>Product</th>
                <th>% Coverage</th>
                <th>% Identity</th>
                <th>Accession</th>
            </tr>
        </thead>
        <tbody>'

    HTML_FOOTER='        </tbody>
    </table>
    <script>
        document.getElementById("analysis-date").textContent = new Date().toLocaleString();
        document.getElementById("database-used").textContent = "~{if use_local_db then "Custom Local Database" else "Standard PlasmidFinder Database"}";
        document.getElementById("samples-processed").textContent = "1";

        // Add coverage highlighting
        document.querySelectorAll("td:nth-child(5)").forEach(cell => {
            const coverage = parseFloat(cell.textContent);
            if (coverage >= 90) cell.classList.add("coverage-high");
            else if (coverage >= 70) cell.classList.add("coverage-medium");
            else cell.classList.add("coverage-low");
        });
    </script>
</body>
</html>'

    # Database setup
    db_to_use="plasmidfinder"
    if [ "~{use_local_db}" == "true" ]; then
      if [ ! -f "~{local_db}" ]; then
        echo "ERROR: Local database file not found" >> mge.log
        echo "MISSING_LOCAL_DB" > mge_results/skipped.txt
        echo "<h1>MGE analysis skipped - local database missing</h1>" > html_results/skipped.html
        exit 0
      fi

      echo "Setting up local MGE database" >> mge.log
      mkdir -p /root/abricate/db/plasmidfinder
      if ! cp "~{local_db}" /root/abricate/db/plasmidfinder/plasmidfinder.fa; then
        echo "ERROR: Failed to copy local database" >> mge.log
        echo "DB_COPY_FAILED" > mge_results/skipped.txt
        echo "<h1>MGE analysis skipped - database setup failed</h1>" > html_results/skipped.html
        exit 0
      fi

      if ! abricate --setupdb --db plasmidfinder --debug >> mge.log 2>&1; then
        echo "ERROR: Failed to setup local MGE database" >> mge.log
        echo "DB_SETUP_FAILED" > mge_results/skipped.txt
        echo "<h1>MGE analysis skipped - database setup failed</h1>" > html_results/skipped.html
        exit 0
      fi

      # Confirmation in log that the run will use the local plasmidfinder
      echo "Database in use: plasmidfinder (LOCAL, from ~{local_db})" >> mge.log
    else
      if ! abricate --list | grep -q "plasmidfinder"; then
        if ! abricate --setupdb --db plasmidfinder >> mge.log 2>&1; then
          echo "ERROR: Failed to setup plasmidfinder database" >> mge.log
          echo "DB_SETUP_FAILED" > mge_results/skipped.txt
          echo "<h1>MGE analysis skipped - database setup failed</h1>" > html_results/skipped.html
          exit 0
        fi
      fi
      echo "Database in use: plasmidfinder (STANDARD)" >> mge.log
    fi

    # Derive sample name and outputs (per-sample only)
    sample_name=$(basename "~{assembly}" | sed -E 's/\.(fa|fasta|fna|fsa|contigs|scaffolds)(\.(gz|bz2))?$//i')
    output_tsv="mge_results/${sample_name}_mge.tsv"
    output_html="html_results/${sample_name}_mge.html"

    echo "Processing ${sample_name}" >> mge.log

    # Run MGE analysis with error handling
    set +e
    abricate \
      --db $db_to_use \
      --mincov 80 \
      --minid 90 \
      --threads ~{cpu} \
      --nopath \
      "~{assembly}" > "${output_tsv}.tmp" 2>> mge.log
    status=$?
    set -e

    if [ $status -ne 0 ] || [ ! -s "${output_tsv}.tmp" ]; then
      echo "WARNING: MGE analysis failed for ${sample_name}" >> mge.log
      echo -e "Sample\tMGE_Type\tGene\tProduct\t%Coverage\t%Identity\tAccession" > "$output_tsv"
      echo -e "${sample_name}\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR" >> "$output_tsv"
    else
      # Process successful output (FIXED AWK: robust header skip + correct column mapping)
      awk -v sample="${sample_name}" '
        BEGIN {
          OFS="\t";
          print "Sample\tMGE_Type\tGene\tProduct\t%Coverage\t%Identity\tAccession";
        }
        NR==1 || $0 ~ /^#/ { next }
        {
          gene=$6;          # correct gene column
          product=$14;      # correct product column
          pcov=$10;         # % coverage
          pid=$11;          # % identity
          accession=$13;    # accession
          print sample, "Plasmid", gene, product, pcov, pid, accession;
        }' "${output_tsv}.tmp" > "$output_tsv"
    fi
    rm -f "${output_tsv}.tmp"

    # Generate HTML report for this sample (template preserved)
    {
      echo "$HTML_HEADER"

      if grep -q "ERROR" "$output_tsv"; then
        echo "<tr class=\"error-row\">"
        echo "<td colspan=\"7\">MGE analysis failed for ${sample_name}</td>"
        echo "</tr>"
      else
        tail -n +2 "$output_tsv" | while IFS=$'\t' read -r _ mge_type gene product coverage identity accession; do
          echo "<tr>"
          echo "<td>${sample_name}</td>"
          echo "<td>${mge_type}</td>"
          echo "<td class=\"mge\">${gene}</td>"
          echo "<td>${product}</td>"
          echo "<td>${coverage}</td>"
          echo "<td>${identity}</td>"
          echo "<td>${accession}</td>"
          echo "</tr>"
        done
      fi

      echo "$HTML_FOOTER"
    } > "$output_html"

    echo "MGE analysis completed at $(date)" >> mge.log
    echo "Output files:" >> mge.log
    ls -lh mge_results/* html_results/* >> mge.log 2>&1 || true
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 75 HDD"
    continueOnReturnCode: true
    timeout: "8 hours"
  }

  output {
    File tsv_out  = if do_mge_analysis then glob("mge_results/*_mge.tsv")[0]  else "mge_results/skipped.txt"
    File html_out = if do_mge_analysis then glob("html_results/*_mge.html")[0] else "html_results/skipped.html"
    File mge_log  = "mge.log"
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
    Boolean skip = false
  }

  command <<<
    set -euo pipefail
    set -x  # Enable debugging

    if [ "~{skip}" == "true" ]; then
      echo "Skipping virulence analysis as requested" > virulence.log
      echo "Creating empty output directories for portability" >> virulence.log

      mkdir -p virulence_results
      mkdir -p html_results

      # No combined outputs; per-sample files intentionally not created.
      echo "Skipping completed at $(date)" >> virulence.log
      exit 0
    fi

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
    echo "- Min coverage: ~{min_coverage}" >> virulence.log
    echo "- Min identity: ~{min_identity}" >> virulence.log
    echo "- CPU: ~{cpu}" >> virulence.log

    mkdir -p virulence_results
    mkdir -p html_results

    # HTML template components
    HTML_HEADER='<!DOCTYPE html>
<html>
<head>
    <title>Virulence Factor Analysis</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; position: sticky; top: 0; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .virulence { font-weight: bold; color: #e74c3c; }
        .coverage-high { background-color: #e6ffe6; }
        .coverage-medium { background-color: #fff9e6; }
        .coverage-low { background-color: #ffe6e6; }
        .risk-high { color: #d63031; }
        .risk-medium { color: #f39c12; }
        .risk-low { color: #27ae60; }
        .summary-card {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
    </style>
</head>
<body>
    <h1>Virulence Factor Analysis Report</h1>
    <div class="summary-card">
        <p><strong>Analysis Date:</strong> <span id="analysis-date"></span></p>
        <p><strong>Database Used:</strong> ~{db_name}</p>
        <p><strong>Parameters:</strong> Min Coverage: ~{min_coverage}%, Min Identity: ~{min_identity}%</p>
    </div>
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Virulence Factor</th>
                <th>Product</th>
                <th>% Coverage</th>
                <th>% Identity</th>
                <th>Risk Level</th>
            </tr>
        </thead>
        <tbody>'

    HTML_FOOTER='        </tbody>
    </table>
    <script>
        document.getElementById("analysis-date").textContent = new Date().toLocaleString();

        // Add coverage highlighting
        const coverageCells = document.querySelectorAll("td:nth-child(4)");
        coverageCells.forEach(cell => {
            const coverage = parseFloat(cell.textContent);
            if (coverage >= 90) cell.classList.add("coverage-high");
            else if (coverage >= 70) cell.classList.add("coverage-medium");
            else cell.classList.add("coverage-low");
        });
    </script>
</body>
</html>'

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

    idx=1
    for asm_file in ~{sep=' ' assembly_output}; do
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//')
      report_prefix="Virulence_Report_${idx}_${sample_name}"
      output_tsv="virulence_results/${report_prefix}.tsv"
      output_html="html_results/${report_prefix}.html"

      echo "Processing $sample_name (Report $idx)" >> virulence.log

      # Run virulence analysis
      abricate \
        --db "$db_to_use" \
        --mincov ~{min_coverage} \
        --minid ~{min_identity} \
        --nopath \
        --quiet \
        --threads ~{cpu} \
        "$asm_file" > "${output_tsv}.tmp" 2>> virulence.log || {
          echo "WARNING: Virulence analysis failed for $sample_name" >> virulence.log
          # Create minimal result file for failed samples
          echo -e "Sample\tVirulence_Factor\tProduct\t%Coverage\t%Identity\tRisk_Level" > "$output_tsv"
          echo -e "$sample_name\tERROR\tERROR\tERROR\tERROR\tERROR" >> "$output_tsv"
          # Still create a minimal HTML to keep downstream stable
          awk -v header="$HTML_HEADER" -v footer="$HTML_FOOTER" -v s="$sample_name" '
          BEGIN { print header }
          {
            print "<tr style='\''background-color: #ffdddd'\''>"
            print "<td>" s "</td><td colspan='\''5'\''>Virulence analysis failed</td>"
            print "</tr>"
          }
          END { print footer }' /dev/null > "$output_html"
          idx=$((idx+1))
          continue
        }

      # Process output to simplified format
      awk -v sample="$sample_name" -v min_cov=~{min_coverage} -v min_id=~{min_identity} '
      BEGIN {
        OFS="\t";
        print "Sample\tVirulence_Factor\tProduct\t%Coverage\t%Identity\tRisk_Level";
      }
      NR==1 { next } # Skip header
      {
        # Abricate columns: 6=GENE, 10=%COVERAGE, 11=%IDENTITY, 14=PRODUCT
        cov = $10 + 0.0;
        idn = $11 + 0.0;

        risk = "Low";
        if (cov >= 90 && idn >= 90)      risk = "High";
        else if (cov >= 70 && idn >= 70) risk = "Medium";

        print sample, $6, $14, cov, idn, risk;
      }' "${output_tsv}.tmp" > "$output_tsv"
      rm "${output_tsv}.tmp"

      # Generate HTML version
      awk -v header="$HTML_HEADER" -v footer="$HTML_FOOTER" -v idx="$idx" -v sample="$sample_name" '
      BEGIN {
        print header;
      }
      NR==1 { next } # Skip header
      {
        risk_class = "risk-low";
        if ($6 == "High") risk_class = "risk-high";
        else if ($6 == "Medium") risk_class = "risk-medium";

        printf "<tr>";
        printf "<td>Virulence Report %d — %s</td>", idx, sample;  # Sample with consistent naming
        printf "<td class=\"virulence\">%s</td>", $2;  # Virulence Factor
        printf "<td>%s</td>", $3;  # Product
        printf "<td>%.1f</td>", $4; # %Coverage
        printf "<td>%.1f</td>", $5; # %Identity
        printf "<td class=\"%s\">%s</td>", risk_class, $6;  # Risk Level
        print "</tr>";
      }
      END {
        print footer;
      }' "$output_tsv" > "$output_html"

      idx=$((idx+1))
    done

    echo "Virulence analysis completed at $(date)" >> virulence.log
    echo "Output files:" >> virulence.log
    ls -lh virulence_results/* html_results/* >> virulence.log 2>&1 || true
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    cpu: cpu
    memory: "6 GB"
    disks: "local-disk 50 HDD"
    preemptible: 2
    continueOnReturnCode: true
    timeout: "8 hours"
  }

  output {
    Array[File] virulence_reports = glob("virulence_results/Virulence_Report_*.tsv")
    Array[File] html_reports      = glob("html_results/Virulence_Report_*.html")
    File        virulence_log     = "virulence.log"
  }
}

task BLAST_ANALYSIS {
  input {
    Array[File]? contig_fastas
    String blast_db = "nt"
    Boolean use_local_blast = false
    File? local_blast_db
    Int max_target_seqs = 250
    Float evalue = 0.000001
    Int min_contig_length = 300
    Int? max_contig_length
    Boolean do_blast = true
    Int cpu = 8
    Int memory_gb = 10
    Boolean skip = false
    Boolean debug = false
    Boolean parse_seqids = false
    Int max_retries_per_sample = 1
    Int retry_delay_seconds = 30
  }

  command <<<
    #!/bin/bash
    set -euo pipefail
    export LC_ALL=C

    # Enhanced debugging
    debug_log() {
      echo "[DEBUG][$(date '+%Y-%m-%d %H:%M:%S')] $1" >> debug.log
    }

    # Initialize execution log
    {
      echo "=== BLAST ANALYSIS STARTED ==="
      echo "Timestamp: $(date +"%Y-%m-%d %H:%M:%S")"
      echo "Input files: ~{sep=' ' contig_fastas}"
      echo "BLAST database: ~{blast_db}"
      echo "Memory allocated: ~{memory_gb}GB"
      echo "Debug mode: ~{debug}"
    } > execution.log

    # Early exit if no contigs provided
    if [ -z "~{sep=' ' contig_fastas}" ]; then
      echo "No contig files provided - creating empty outputs" >> execution.log
      mkdir -p blast_results
      touch sample_ids.txt
      exit 0
    fi

    # Create directories
    mkdir -p blast_results parallel_tmp || {
      echo "ERROR: Failed to create output directories" >> execution.log
      exit 1
    }
    export TMPDIR
    TMPDIR=$(mktemp -d ./blast_tmp.XXXXXX)
    debug_log "Created temp directory: $TMPDIR"

    # Handle BLAST database setup
    if [ "~{use_local_blast}" = "true" ] && [ -n "~{local_blast_db}" ]; then
      echo "Setting up local BLAST database from input file" >> execution.log
      makeblastdb -in "~{local_blast_db}" -dbtype nucl -parse_seqids \
                  -title "~{blast_db}" -out "~{blast_db}" \
                  -logfile "makeblastdb.log" || {
        echo "Failed to build local DB" >> execution.log
        exit 1
      }
      export BLAST_DB="~{blast_db}"
    else
      echo "Using remote or prebuilt DB: ~{blast_db}" >> execution.log
      export BLAST_DB="~{blast_db}"
    fi

    # System diagnostics
    {
      echo -e "\n=== SYSTEM CHECK ==="
      echo "CPU cores: $(nproc)"
      echo "Memory: $(free -h | awk '/Mem:/{print $2}') available"
      echo "Disk space: $(df -h . | awk 'NR==2{print $4}') free"
      echo "BLAST version: $(blastn -version 2>&1 | head -1)"
      echo "Python version: $(python3 --version 2>&1)"
      if command -v parallel >/dev/null 2>&1; then
        echo "GNU parallel: $(parallel --version | head -1)"
      else
        echo "WARNING: GNU parallel not found in PATH"
      fi
    } >> execution.log

    # HTML generator: reads top_hits.tsv ONLY (no species-summary block)
    cat > tsv_to_html.py <<'PYEOF'
import sys
import os
import re
from datetime import datetime

CSS_STYLE = """
<style>
  body { font-family: Arial, Helvetica, sans-serif; margin: 20px; color: #333; line-height: 1.6; }
  .container { max-width: 1200px; margin: 0 auto; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  table { border-collapse: collapse; width: 100%; margin: 20px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
  th { background-color: #3498db; color: white; text-align: left; padding: 12px; position: sticky; top: 0; }
  td { padding: 10px; border-bottom: 1px solid #ddd; }
  tr:nth-child(even) { background-color: #f9f9f9; }
  tr:hover { background-color: #eaf2f8; }
  .footer { margin-top: 30px; font-size: 0.9em; color: #7f8c8d; text-align: right; }
  .no-results { color: #e74c3c; font-style: italic; }
  .error { background-color: #fdecea; padding: 15px; border-radius: 4px; border-left: 4px solid #e74c3c; }
</style>
"""

HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>BLAST Report - {sample_id}</title>
  {css}
</head>
<body>
  <div class="container">
    <h1>BLAST Results: {sample_id}</h1>

    {content}

    <div class="footer">
      <p>Report generated: {timestamp}</p>
      <p>Database: {db_name}</p>
    </div>
  </div>
</body>
</html>"""

# Bold+italicize the Latin binomial (Genus species) at the start of the species string
_BINOMIAL_RE = re.compile(r"^([A-Z][a-z]+ [a-z]+)(.*)$")
def _format_species_cell(text: str) -> str:
  if not text:
    return text
  m = _BINOMIAL_RE.match(text)
  if m:
    return f"<i><b>{m.group(1)}</b></i>{m.group(2)}"
  return text

def generate_report(sample_dir):
  sample_id = os.path.basename(sample_dir)
  db_name = os.environ.get("BLAST_DB", "nt")
  input_tsv = os.path.join(sample_dir, "top_hits.tsv")
  output_html = os.path.join(sample_dir, f"BLAST_Report_{sample_id}.html")
  timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

  if (not os.path.exists(input_tsv)) or os.path.getsize(input_tsv) <= 1:
    content = '<p class="no-results">No significant BLAST hits found</p>'
  else:
    try:
      with open(input_tsv, "r", encoding="utf-8") as f:
        headers = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f]

      # Display header: rename stitle -> species (safety even if upstream already renamed)
      headers = ["species" if h == "stitle" else h for h in headers]

      # Find the species column index (if present)
      species_idx = None
      for i, h in enumerate(headers):
        if h == "species":
          species_idx = i
          break

      # Build HTML table
      parts = []
      parts.append("<table>")
      parts.append("  <thead><tr>")
      parts.extend(f"    <th>{h}</th>" for h in headers)
      parts.append("  </tr></thead>")
      parts.append("  <tbody>")

      for row in rows:
        parts.append("    <tr>")
        for i, cell in enumerate(row):
          if species_idx is not None and i == species_idx:
            cell = _format_species_cell(cell)
          parts.append(f"      <td>{cell}</td>")
        parts.append("    </tr>")

      parts.append("  </tbody>")
      parts.append("</table>")
      content = "\n".join(parts)

    except Exception as e:
      content = """
      <div class="error">
        <h2>Processing Error</h2>
        <p>Failed to generate report:</p>
        <pre>{}</pre>
      </div>
      """.format(str(e))

  with open(output_html, "w", encoding="utf-8") as f:
    f.write(HTML_TEMPLATE.format(
      sample_id=sample_id,
      css=CSS_STYLE,
      content=content,
      timestamp=timestamp,
      db_name=db_name
    ))

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: python3 tsv_to_html.py <sample_directory>")
    sys.exit(1)
  generate_report(sys.argv[1])
PYEOF

    # Prepare safe working copies of contig FASTAs (avoid writing to read-only inputs)
    work_list="work_list.txt"
    : > "$work_list"
    idx=1
    for fa in ~{sep=' ' contig_fastas}; do
      bn="$(basename "$fa")"
      wf="$TMPDIR/${bn}"
      cp -f -- "$fa" "$wf"
      chmod u+rw -- "$wf" || true
      echo "$wf" >> "$work_list"
      idx=$((idx+1))
    done

    # Main processing function
    process_sample() {
      local contig_file="$1"
      local sample_id
      sample_id=$(basename "$contig_file" | cut -d'.' -f1)
      local sample_dir="blast_results/${sample_id}"

      mkdir -p "$sample_dir"
      debug_log "Processing sample: $sample_id"

      # 1. Filter contigs by length
      awk -v min_len=~{min_contig_length} \
      -v max_len="~{if defined(max_contig_length) then max_contig_length else "NULL"}" \
          'BEGIN {RS=">"; FS="\n"}
          NR>1 {
            seq="";
            for(i=2;i<=NF;i++) seq=seq $i;
            if(length(seq)>=min_len && (max_len == "NULL" || length(seq)<=max_len)) {
              print ">"$1;
              for(i=2;i<=NF;i++) print $i
            }
          }' "$contig_file" > "${sample_dir}/filtered.fa"

      debug_log "Filtered contigs: $(grep -c '^>' "${sample_dir}/filtered.fa" || true) sequences"

      # 2. Run BLAST if we have sequences
      if [ -s "${sample_dir}/filtered.fa" ]; then
        blastn \
          -query "${sample_dir}/filtered.fa" \
          -db "${BLAST_DB}" \
          -outfmt "6 std qlen slen stitle" \
          -out "${sample_dir}/blast_results.tsv" \
          -evalue ~{evalue} \
          -max_target_seqs ~{max_target_seqs} \
          -num_threads 1 \
          -task megablast > "${sample_dir}/blast.log" 2>&1 || true
      else
        echo "No sequences passed length filters" > "${sample_dir}/blast.log"
        : > "${sample_dir}/blast_results.tsv"
      fi

      debug_log "BLAST results: $(wc -l < "${sample_dir}/blast_results.tsv") hits"

      # 3. Select top hits (by bitscore), and write header with 'species' instead of 'stitle'
      {
        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tspecies"
        if [ -s "${sample_dir}/blast_results.tsv" ]; then
          sort -t$'\t' -k12,12gr "${sample_dir}/blast_results.tsv" | head -10 | \
          awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}'
        fi
      } > "${sample_dir}/top_hits.tsv"

      debug_log "Top hits (renamed stitle->species): $(wc -l < "${sample_dir}/top_hits.tsv") records"

      # 4. Generate HTML report (NO species-summary block, species is bold+italic)
      python3 tsv_to_html.py "${sample_dir}" > "${sample_dir}/html_generation.log" 2>&1 || {
        echo "<html><body><h1>Report Generation Error</h1><pre>$(cat ${sample_dir}/html_generation.log)</pre></body></html>" \
          > "${sample_dir}/BLAST_Report_${sample_id}.html"
      }

      echo "$sample_id"
    }
    export -f process_sample debug_log
    export BLAST_DB

    # Process samples in parallel from the working-copy list
    /usr/bin/parallel --will-cite --noswap -j ~{cpu} --joblog parallel.log \
      --progress --eta --tagstring "{}" "process_sample {}" :::: "$work_list" \
      1> sample_ids.txt


    # Finalization with results summary
    {
      echo -e "\n=== PROCESSING COMPLETED ==="
      echo "Successful samples: $(grep -c . sample_ids.txt || echo 0)"
      echo "BLAST database used: $BLAST_DB"
      echo "Final output checks:"
      find blast_results -name "top_hits.tsv" -exec sh -c 'echo "  {}: $(($(wc -l < "{}")-1)) records"' \;
      echo "Timestamp: $(date +"%Y-%m-%d %H:%M:%S")"
    } >> execution.log

    # Cleanup
    rm -rf "${TMPDIR}"
  >>>
  runtime {
    docker: "gmboowa/blast-analysis:1.9.4"
    cpu: select_first([cpu, 1])
    memory: "~{memory_gb} GB"
    disks: "local-disk 100 HDD"
    preemptible: 2
    dockerOptions: "--entrypoint /bin/bash --user 0:0"
  }

  output {
    Array[File] blast_top10 = glob("blast_results/*/top_hits.tsv")
    Array[File] blast_logs = glob("blast_results/*/blast.log")
    Array[File] blast_reports = glob("blast_results/*/BLAST_Report_*.html")
    Array[File] blast_results = glob("blast_results/*/blast_results.tsv")
    Array[File] filtered_contigs = glob("blast_results/*/filtered.fa")
    Array[String] sample_ids = read_lines("sample_ids.txt")
    File execution_log = "execution.log"
    File parallel_log = "parallel.log"
    File? makeblastdb_log = "makeblastdb.log"
    File? debug_log = "debug.log"
  }
}
task TREE_VISUALIZATION {
  input {
    File input_tree
    Int width = 1600
    Int height = 1200
    String image_format = "png"
    Boolean force_offscreen = true
    String? tree_title
    Int title_font_size = 18
    String layout = "rectangular"
    Boolean show_branch_lengths = false
    Boolean show_scale = true
    Float branch_thickness = 3.0
    String color_scheme = "standard"
    Int label_font_size = 14
    Int? font_size
    Boolean label_bold = true
    String label_color = "#111111"
    Boolean label_heavy = false
    Int tip_size = 9
    String tip_fill_color = "#00008B"
    String tip_border_color = "#00008B"
    String background_color = "#FFFFFF"
    Int margins_px = 48
    String line_cap = "round"
    Boolean autoscale_with_thickness = true
    Boolean bootstrap_tools = false
  }

  command <<<
    set -euo pipefail
    : > debug.log
    echo "[INFO] TREE_VISUALIZATION starting $(date)" | tee -a debug.log
    if [ "~{force_offscreen}" = "true" ]; then
      export QT_QPA_PLATFORM=offscreen
      export MPLBACKEND=Agg
      echo "[INFO] Offscreen rendering enabled" | tee -a debug.log
    fi
    mkdir -p output
    TREE_BASE="$(basename "~{input_tree}" .nwk)"
    OUT_IMG="output/${TREE_BASE}.~{image_format}"
    OUT_LOG="output/${TREE_BASE}.log"
    ERR_LOG="output/${TREE_BASE}_error.log"
    echo "[INFO] Input tree: ~{input_tree} ($(wc -c < "~{input_tree}") bytes)" | tee -a debug.log
    if [ "~{bootstrap_tools}" = "true" ]; then
      echo "[BOOTSTRAP] starting…" | tee -a debug.log
      mkdir -p .mamba
      export MAMBA_ROOT_PREFIX="$PWD/.mamba"
      export MAMBA_EXE="$PWD/.mamba/micromamba"
      if [ ! -x "$MAMBA_EXE" ]; then
        echo "[BOOTSTRAP] fetching micromamba…" | tee -a debug.log
        (curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C .mamba --strip-components=1 bin/micromamba) \
        || (wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C .mamba --strip-components=1 bin/micromamba) \
        || echo "[BOOTSTRAP] micromamba download failed" | tee -a debug.log
      fi
      if [ -x "$MAMBA_EXE" ]; then
        set +e
        "$MAMBA_EXE" create -y -n viz -c conda-forge \
          python=3.11 ete3 pyqt pillow matplotlib biopython \
          cairo pango pixman freetype fontconfig harfbuzz cairocffi \
          > .mamba/create.log 2>&1
        rc=$?
        set -e
        if [ $rc -eq 0 ]; then
          export PATH="$PWD/.mamba/envs/viz/bin:$PATH"
          echo "[BOOTSTRAP] env ready" | tee -a debug.log
        else
          echo "[BOOTSTRAP] env create failed (rc=$rc); continuing with base image" | tee -a debug.log
        fi
      fi
    fi
    if [ $(wc -c < "~{input_tree}") -le 10 ]; then
      echo "[ERROR] Input tree appears empty; writing placeholder" | tee -a debug.log
      python3 - <<'PY'
from PIL import Image, ImageDraw
w,h=~{width},~{height}
img=Image.new('RGB',(w,h),'white')
ImageDraw.Draw(img).text((20,20),"Empty tree input",fill='black')
img.save("output/~{basename(input_tree, '.nwk')}.~{image_format}")
PY
      echo "Empty/invalid tree" > "${ERR_LOG}"
      echo "Rendered placeholder at $(date)" > "${OUT_LOG}"
      exit 0
    fi
    set +e
    python3 - <<'PY'
import colorsys, random
from datetime import datetime
try:
    from ete3 import Tree, TreeStyle, TextFace
except Exception as e:
    open("output/~{basename(input_tree, '.nwk')}_error.log","a").write(f"[ete3 import error] {e}\n")
    raise SystemExit(20)
tree_path  = "~{input_tree}"
title      = "~{tree_title}" if "~{tree_title}" != "None" else ""
layout     = "~{layout}".lower()
show_bl    = ~{if show_branch_lengths then "True" else "False"}
show_scale = ~{if show_scale then "True" else "False"}
bg         = "~{background_color}"
margins    = ~{margins_px}
thick      = float(~{branch_thickness})
color_mode = "~{color_scheme}".lower()
base_label_size = int(~{select_first([font_size, label_font_size])})
label_bold      = ~{if label_bold then "True" else "False"}
label_color     = "~{label_color}"
label_heavy     = ~{if label_heavy then "True" else "False"}
base_tip_size   = int(~{tip_size})
tip_fill        = "~{tip_fill_color}"
tip_border      = "~{tip_border_color}"
autoscale       = ~{if autoscale_with_thickness then "True" else "False"}
if autoscale:
    tip_size = max(5, min(22, int(round(thick * 2.5))))
    label_size = max(12, min(28, int(round(10 + thick * 1.6))))
else:
    tip_size = base_tip_size
    label_size = base_label_size
try:
    t = Tree(tree_path)
except Exception as e:
    open("output/~{basename(input_tree, '.nwk')}_error.log","a").write(f"[ete3 load error] {e}\n")
    raise SystemExit(21)
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = show_bl
ts.show_scale = show_scale
ts.branch_vertical_margin = 14
ts.scale = max(1.0, ~{width} / 6.0)
ts.margin_left = margins
ts.margin_right = margins
ts.margin_top = margins
ts.margin_bottom = margins
ts.bgcolor = bg
if title:
    ts.title.add_face(TextFace(title, fsize=~{title_font_size}, bold=True), column=0)
if layout == "circular":
    ts.mode = "c"; ts.arc_start = 0; ts.arc_span = 360
elif layout == "fan":
    ts.mode = "c"; ts.arc_start = -180; ts.arc_span = 240
else:
    ts.mode = "r"; ts.root_opening_factor = 1.0
if color_mode == "gradient":
    max_dist = max([n.dist for n in t.traverse()], default=1.0) or 1.0
    def color_for(n):
        ratio = 0.0 if max_dist==0 else min(max(n.dist/max_dist,0.0),1.0)
        r,g,b = colorsys.hsv_to_rgb(0.66*(1-ratio), 0.65, 0.92)
        return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
elif color_mode == "categorical":
    random.seed(42)
    hues = {}
    def color_for(n):
        parent = n.up
        if parent not in hues:
            hues[parent] = random.uniform(0.0, 0.8)
        r,g,b = colorsys.hsv_to_rgb(hues[parent], 0.8, 0.90)
        return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
else:
    def color_for(n): return "#000000"
for n in t.traverse():
    col = color_for(n)
    n.img_style["hz_line_width"] = thick
    n.img_style["vt_line_width"] = thick
    n.img_style["hz_line_color"] = col
    n.img_style["vt_line_color"] = col
    if n.is_leaf():
        n.img_style["size"] = max(1, tip_size)
        n.img_style["shape"] = "circle"
        n.img_style["fgcolor"] = tip_border
        n.img_style["bgcolor"] = tip_fill
        tf_main = TextFace(n.name, fsize=label_size, fgcolor=label_color, bold=label_bold)
        if label_heavy:
            tf_shadow = TextFace(n.name, fsize=label_size, fgcolor=label_color, bold=True)
            n.add_face(tf_shadow, column=0, position="aligned")
        n.add_face(tf_main, column=0, position="aligned")
    else:
        n.img_style["size"] = 0
out_img = "output/~{basename(input_tree, '.nwk')}.~{image_format}"
try:
    t.render(out_img, w=~{width}, h=~{height}, units="px", tree_style=ts, dpi=300)
    open("output/~{basename(input_tree, '.nwk')}.log","w").write(f"ete3 Rendered OK at {datetime.now().isoformat()}\n")
    raise SystemExit(0)
except Exception as e:
    open("output/~{basename(input_tree, '.nwk')}_error.log","a").write(f"[ete3 render error] {e}\n")
    raise SystemExit(22)
PY
    ETE_RC=$?
    set -e
    echo "[INFO] ete3 render exit code: ${ETE_RC}" | tee -a debug.log
    if [ $ETE_RC -ne 0 ]; then
      echo "[WARN] Falling back to Bio.Phylo + matplotlib" | tee -a debug.log
      set +e
      python3 - <<'PY'
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
tree_path="~{input_tree}"
out_img="output/~{basename(input_tree, '.nwk')}.~{image_format}"
w,h=~{width},~{height}
dpi=150.0
fig = plt.figure(figsize=(w/dpi, h/dpi), dpi=dpi, facecolor="~{background_color}")
ax = fig.add_subplot(1,1,1, facecolor="~{background_color}")
tree = Phylo.read(tree_path,"newick")
Phylo.draw(tree, do_show=False, axes=ax)
for line in ax.get_lines():
    try:
        line.set_linewidth(~{branch_thickness})
        line.set_solid_capstyle("~{line_cap}")
    except Exception:
        pass
for text in ax.texts:
    text.set_fontsize(~{select_first([font_size, label_font_size])})
    text.set_color("~{label_color}")
    text.set_fontweight("bold" if ~{if label_bold then "True" else "False"} else "normal")
if "~{tree_title}" != "None" and len("~{tree_title}")>0:
    ax.set_title("~{tree_title}", fontsize=~{title_font_size}, fontweight='bold', color="#111111")
ax.set_axis_off()
plt.margins(x=0.02, y=0.02)
plt.tight_layout(pad=~{margins_px}/100)
fig.savefig(out_img, dpi=dpi, facecolor="~{background_color}")
open("output/~{basename(input_tree, '.nwk')}.log","w").write(f"Phylo Rendered OK at {datetime.now().isoformat()}\n")
PY
      PHYLO_RC=$?
      set -e
      echo "[INFO] Phylo/matplotlib render exit code: ${PHYLO_RC}" | tee -a debug.log
      if [ $PHYLO_RC -ne 0 ]; then
        echo "[ERROR] All renderers failed; no image produced" | tee -a debug.log
        printf "All renderers failed (ete3=%s, phylo=%s)\n" "$ETE_RC" "$PHYLO_RC" > "${ERR_LOG}"
        exit 0
      fi
    fi
    echo "[INFO] TREE_VISUALIZATION finished $(date)" | tee -a debug.log
  >>>

  runtime {
    docker: "gmboowa/ete3-render:1.18"
    memory: "16 GB"
    cpu: 2
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
    preemptible: 2
  }

  output {
    File final_image = "output/~{basename(input_tree, '.nwk')}.~{image_format}"
    File render_log  = "output/~{basename(input_tree, '.nwk')}.log"
    File debug_log   = "debug.log"
    File? error_log  = "output/~{basename(input_tree, '.nwk')}_error.log"
    Boolean tree_valid  = size(input_tree, "B") > 10
  }
}


task MERGE_REPORTS {
  input {
    # Section HTMLs (optional)
    File? trimming_report_html
    File? assembly_stats_html
    File? annotation_summary_html
    File? pangenome_report_html
    File? gene_heatmap_png
    File? mlst_combined_html
    File? variant_summary_html

    # Per-sample arrays
    Array[File] amr_reports = []
    Array[File] mge_reports = []
    Array[File] virulence_reports = []
    Array[File] quality_reports = []
    Array[File] blast_reports = []

    # Direct-file inputs
    File? pangenome_report
    File? gene_heatmap

    # Collections
    Array[File] tree_images = []

    # Meta
    String workflow_name = "rMAP Analysis Pipeline"
    String version = "2.0"
    String footer_sentence = ""
    Boolean skip = false

    # Logo file as input (optional override)
    File? logo_file
  }

  command <<<
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob

    LOG_FILE="report_generation.log"
    exec > >(tee -a "$LOG_FILE") 2>&1
    echo "=== Starting report generation at $(date) ==="

    : > skip_report.log

    echo "Cleaning working directory..."
    rm -rfv final_report default_report.tgz final_report.tgz || true
    mkdir -pv final_report/assets/images final_report/assets/sections

    readonly REPORT_FILENAME="rMAP_~{version}_final_report.html"
    readonly RUN_DATE="$(date +"%Y-%m-%d %H:%M:%S %Z")"
    export RUN_DATE  # ensure available to tools like perl/sed if needed

    LOGO_HTML_TOP=""
    LOGO_HTML_BOTTOM=""
    LOGO_PATH="final_report/assets/images/rMAP_logo.png"

    # ----------------------------------------------------------
    # Robust logo setup: copy if provided, otherwise download it
    # Try 5 different methods in order until one succeeds.
    # ----------------------------------------------------------
    copy_provided_logo() {
      if [[ -n "~{logo_file}" && -f "~{logo_file}" ]]; then
        echo "[logo] Using provided logo file: ~{logo_file}"
        cp -v "~{logo_file}" "$LOGO_PATH" && return 0 || true
      fi
      return 1
    }

    ensure_tools() {
      # Install minimal tools if missing
      if ! command -v curl >/dev/null 2>&1 || ! command -v wget >/dev/null 2>&1 || ! command -v git >/dev/null 2>&1 || ! command -v python3 >/dev/null 2>&1; then
        echo "[deps] Installing missing tools (curl, wget, git, python3, ca-certificates)..."
        (apt-get update && apt-get install -y --no-install-recommends curl wget git python3 ca-certificates) || echo "[deps] Warning: apt-get failed; continuing with available tools" >&2
      fi
    }

    try_curl_raw() {
      local url="https://raw.githubusercontent.com/gmboowa/rMAP-2.0/main/rMAP-2.0_logo.png"
      if command -v curl >/dev/null 2>&1; then
        echo "[logo] curl raw.githubusercontent attempt..."
        curl -fSL "$url" -o "$LOGO_PATH" && return 0 || true
      fi
      return 1
    }

    try_wget_raw() {
      local url="https://raw.githubusercontent.com/gmboowa/rMAP-2.0/main/rMAP-2.0_logo.png"
      if command -v wget >/dev/null 2>&1; then
        echo "[logo] wget raw.githubusercontent attempt..."
        wget -qO "$LOGO_PATH" "$url" && return 0 || true
      fi
      return 1
    }

    try_github_api() {
      # Use GitHub API with Accept header for raw content
      local api_url="https://api.github.com/repos/gmboowa/rMAP-2.0/contents/rMAP-2.0_logo.png"
      if command -v curl >/dev/null 2>&1; then
        echo "[logo] curl GitHub API attempt..."
        curl -fSL -H "Accept: application/vnd.github.v3.raw" "$api_url" -o "$LOGO_PATH" && return 0 || true
      fi
      return 1
    }

    try_jsdelivr() {
      # CDN mirror via jsDelivr
      local cdn_url="https://cdn.jsdelivr.net/gh/gmboowa/rMAP-2.0@main/rMAP-2.0_logo.png"
      if command -v curl >/dev/null 2>&1; then
        echo "[logo] curl jsDelivr attempt..."
        curl -fSL "$cdn_url" -o "$LOGO_PATH" && return 0 || true
      fi
      if command -v wget >/dev/null 2>&1; then
        echo "[logo] wget jsDelivr attempt..."
        wget -qO "$LOGO_PATH" "$cdn_url" && return 0 || true
      fi
      return 1
    }

    try_git_clone() {
      # Shallow clone and copy the file
      if command -v git >/dev/null 2>&1; then
        echo "[logo] git clone attempt..."
        rm -rf /tmp/rmap2logo || true
        git clone --depth 1 https://github.com/gmboowa/rMAP-2.0.git /tmp/rmap2logo && \
          cp -v /tmp/rmap2logo/rMAP-2.0_logo.png "$LOGO_PATH" && { rm -rf /tmp/rmap2logo || true; return 0; }
        rm -rf /tmp/rmap2logo || true
      fi
      return 1
    }

    try_python_urllib() {
      if command -v python3 >/dev/null 2>&1; then
        echo "[logo] python urllib attempt..."
        python3 - <<'PY' || exit 1
import sys, ssl, urllib.request
url = "https://raw.githubusercontent.com/gmboowa/rMAP-2.0/main/rMAP-2.0_logo.png"
ctx = ssl.create_default_context()
try:
    with urllib.request.urlopen(url, context=ctx, timeout=60) as r, open("final_report/assets/images/rMAP_logo.png","wb") as f:
        f.write(r.read())
    sys.exit(0)
except Exception as e:
    sys.exit(1)
PY
        return 0
      fi
      return 1
    }

    logo_ok=false
    copy_provided_logo && logo_ok=true || true
    if [[ "$logo_ok" != true ]]; then ensure_tools || true; fi
    if [[ "$logo_ok" != true ]]; then try_curl_raw   && logo_ok=true || true; fi
    if [[ "$logo_ok" != true ]]; then try_wget_raw   && logo_ok=true || true; fi
    if [[ "$logo_ok" != true ]]; then try_github_api && logo_ok=true || true; fi
    if [[ "$logo_ok" != true ]]; then try_jsdelivr   && logo_ok=true || true; fi
    if [[ "$logo_ok" != true ]]; then try_git_clone  && logo_ok=true || true; fi
    if [[ "$logo_ok" != true ]]; then try_python_urllib && logo_ok=true || true; fi

    if [[ "$logo_ok" == true && -s "$LOGO_PATH" ]]; then
      echo "[logo] Logo obtained successfully: $LOGO_PATH"
    else
      echo "[logo] WARNING: Failed to obtain logo after multiple methods; captions will render, image will hide." >&2
    fi

    # --- Always build the logo snippet (image hides itself if missing) ---
    cat <<EOT > "/tmp/logo_snippet.html"
<div class="logo-combo">
  <div class="logo-caption-top">rMAP-2.0 Comprehensive Analysis Report Summary: ${RUN_DATE}</div>
  <div class="logo-container"><img src="assets/images/rMAP_logo.png" alt="rMAP Logo" style="height:200px;object-fit:contain;margin:10px" onerror="this.style.display='none'"></div>
  <div class="logo-caption-bottom">This run was generated using rMAP version ~{version} performed on ${RUN_DATE}... Thank you for using this pipeline!!</div>
</div>
EOT
    LOGO_HTML_TOP="$(cat /tmp/logo_snippet.html)"
    LOGO_HTML_BOTTOM="$LOGO_HTML_TOP"

    if [[ "~{skip}" == "true" ]]; then
      echo "Skipping report generation as requested" > skip_report.log
      cat > "final_report/${REPORT_FILENAME}" <<EOF
<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Report Generation Skipped</title>
<style>
  body { font-family: Arial, sans-serif; margin: 40px; }
  .skip-notice { background: #fff3cd; border: 1px solid #ffeeba; padding: 20px; border-radius: 5px; margin: 20px 0; text-align: center; }
  .logo-container { text-align: center; margin: 10px 0; }
  .logo-caption-top, .logo-caption-bottom { text-align: center; margin: 6px 0; font-weight: 600; }
</style></head>
<body>
  ${LOGO_HTML_TOP}
  <h1>Report Generation Skipped</h1>
  <div class="skip-notice">
    <p>This report was not generated because the workflow was configured to skip this step.</p>
  </div>
  ${LOGO_HTML_BOTTOM}
</body></html>
EOF
      tar -czvf final_report.tgz final_report
      exit 0
    fi

    echo "Creating default report (fallback)..."
    cat > "final_report/${REPORT_FILENAME}" <<EOF
<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Report Generation Failed</title>
<style>.logo-container{text-align:center;margin:10px 0}.logo-caption-top,.logo-caption-bottom{text-align:center;margin:6px 0;font-weight:600}</style>
</head><body>
  ${LOGO_HTML_TOP}
  <h1>Report Generation Failed</h1>
  <p>An error occurred while generating the report. Please check the workflow logs.</p>
  ${LOGO_HTML_BOTTOM}
</body></html>
EOF
    tar -czvf default_report.tgz "final_report/${REPORT_FILENAME}"

    safe_copy() {
      local src="$1"; local dst="$2"
      [[ -n "${src}" && -f "${src}" ]] || { echo "Skip copy: missing ${src}" >&2; return 1; }
      mkdir -p "$(dirname "${dst}")"
      cp -v "${src}" "${dst}" || { echo "Copy failed: ${src}" >&2; return 1; }
    }

    # ----------------------------- ROBUST sample ID extraction -----------------------------
    extract_sample_id() {
      local path="$1"
      local fname="$(basename "$path")"
      local stem="${fname%.*}"

      # 1) Prefer an accession anywhere in the full path
      local acc
      acc="$(printf '%s\n' "$path" | grep -Eo '(SRR|ERR|DRR)[0-9]+' | head -n1 || true)"
      if [[ -n "${acc:-}" ]]; then
        echo "$acc"
        return 0
      fi

      # 2) Prefer tokens like A55870 / KPN123 etc. (letters+digits)
      local alnum
      alnum="$(printf '%s\n' "$path" | grep -Eo '[A-Za-z][A-Za-z0-9]*[0-9]{2,}' | tail -n1 || true)"
      if [[ -n "${alnum:-}" && ! "${alnum,,}" =~ ^(report|blast|virulence|amr|mge)$ ]]; then
        echo "$alnum"
        return 0
      fi

      # 3) Walk up directories and filename stem; clean and try again
      local d1="$(basename "$(dirname "$path")")"
      local d2="$(basename "$(dirname "$(dirname "$path")")")"
      local d3="$(basename "$(dirname "$(dirname "$(dirname "$path")")")")"

      is_bad_token() {
        local t="$1"
        [[ -z "$t" || "$t" =~ ^-?[0-9]+$ || "$t" =~ ^(execution|inputs|attempt-[0-9]+|shard-[0-9]+|glob-.*|tmp.*|call-[^/]+|combined)$ ]]
      }

      clean_token() {
        echo "$1" | sed -E '
          s/\.(html?|txt)$//;
          s/(^|[_-])(amr|mge|virulence|blast|report)([_-]|$)/_/Ig;
          s/^(shard-|scatter-|attempt-|call-)?[0-9]+[_-]+//I;
          s/[_-]+$//;
        '
      }

      local cand acc2
      for raw in "$d1" "$d2" "$d3" "$stem"; do
        cand="$(clean_token "$raw")"
        acc2="$(printf '%s\n' "$cand" | grep -Eo '(SRR|ERR|DRR)[0-9]+' | head -n1 || true)"
        if [[ -n "${acc2:-}" ]]; then
          echo "$acc2"; return 0
        fi
        alnum="$(printf '%s\n' "$cand" | grep -Eo '[A-Za-z][A-Za-z0-9]*[0-9]{2,}' | tail -n1 || true)"
        if [[ -n "${alnum:-}" && ! "${alnum,,}" =~ ^(report|blast|virulence|amr|mge)$ ]]; then
          echo "$alnum"; return 0
        fi
        if ! is_bad_token "$cand"; then
          echo "$cand"; return 0
        fi
      done

      echo "${stem:-sample}"
    }
    # --------------------------------------------------------------------------------------

    embed_html() {
      local label="$1"; local file="$2"; local out="$3"
      {
        echo "          <div class=\"report-card\">"
        echo "            <h3>${label}</h3>"
        echo "            <div class=\"report-content\">"
        if [[ -f "${file}" ]]; then
          if grep -q "<table" "${file}" 2>/dev/null; then
            sed -n '/<table[ >]/,/<\/table>/p' "${file}"
          else
            sed 's/<script[^>]*>.*<\/script>//gI; /<link[^>]*>/d' "${file}"
          fi
        else
          echo "<p><em>No ${label} data available.</em></p>"
        fi
        echo "            </div>"
        echo "          </div>"
      } >> "${out}"
    }

    # --- Modified Heatmap Handling with Explicit Precedence ---
    GH_BASENAME=""
    if [[ -n "~{gene_heatmap}" && -f "~{gene_heatmap}" ]]; then
      echo "[heatmap] Using primary gene_heatmap input"
      if safe_copy "~{gene_heatmap}" "final_report/assets/images/$(basename "~{gene_heatmap}")"; then
        GH_BASENAME="$(basename "~{gene_heatmap}")"
      fi
    elif [[ -n "~{gene_heatmap_png}" && -f "~{gene_heatmap_png}" ]]; then
      echo "[heatmap] Falling back to gene_heatmap_png input"
      if safe_copy "~{gene_heatmap_png}" "final_report/assets/images/$(basename "~{gene_heatmap_png}")"; then
        GH_BASENAME="$(basename "~{gene_heatmap_png}")"
      fi
    fi

    if [[ -n "~{gene_heatmap}" && -f "~{gene_heatmap}" && -n "~{gene_heatmap_png}" && -f "~{gene_heatmap_png}" ]]; then
      echo "WARNING: Both gene_heatmap and gene_heatmap_png provided - using gene_heatmap" >&2
    fi

    # Copy tree images up front (if any)
    for img in ~{sep=' ' tree_images}; do
      safe_copy "${img}" "final_report/assets/images/$(basename "${img}")" || true
    done

    echo "Generating main HTML report..."
    cat > "final_report/${REPORT_FILENAME}" <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>~{workflow_name} - Comprehensive Report</title>
  <style>
    :root { --primary:#3498db; --secondary:#2c3e50; --light:#f8f9fa; --text:#333; --accent:#fff3cd; }
    body { font-family:'Segoe UI',Tahoma,Arial,sans-serif; margin:0; background:#f5f5f5; color:var(--text); line-height:1.6; }
    .container { max-width:1200px; margin:0 auto; padding:20px; }
    .header { background:var(--secondary); color:#fff; padding:24px; border-radius:6px; margin-bottom:24px; }
    .logo-combo { margin: 20px 0; }
    .logo-container { text-align:center; margin:10px 0 }
    .logo-caption-top, .logo-caption-bottom { text-align:center; margin:6px 0; font-weight:600 }
    h1,h2,h3{margin-top:0} h1{font-size:2rem;text-align:center}
    h2{font-size:1.75rem;color:var(--secondary); border-bottom:2px solid #3498db; padding-bottom:8px;}
    h3{font-size:1.25rem;color:var(--secondary)}
    .meta{opacity:.9;text-align:center;font-size:.9rem}
    .layout{display:flex;gap:24px}
    .nav{position:sticky;top:16px;width:260px;background:#f8f9fa;padding:16px;border-radius:6px;height:fit-content}
    .nav a{display:block;padding:8px 10px;border-radius:4px;text-decoration:none;color:#2c3e50;margin:4px 0}
    .nav a:hover{background:#3498db;color:#fff}
    .content{flex:1;background:#fff;padding:20px;border-radius:6px;box-shadow:0 2px 4px rgba(0,0,0,.05)}
    .section{margin-bottom:36px;border-bottom:1px solid #eee;padding-bottom:20px}.section:last-child{border:none}
    .report-card{background:#f8f9fa;padding:16px;border-radius:6px;margin:12px 0}
    .image-center{text-align:center;margin:12px 0}
    img{max-width:100%;height:auto;border:1px solid #ddd;border-radius:4px}
    table{width:100%;border-collapse:collapse;margin:12px 0}
    th,td{padding:8px 12px;text-align:left;border:1px solid #ddd}
    th{background-color:#2c3e50;color:white} tr:nth-child(even){background:#f2f2f2}
    .footer-banner{text-align:center;color:#333;margin:24px 0;padding:10px 12px;background:#fff3cd;border:1px solid #e6cf8b;border-radius:6px;font-weight:600}
    @media(max-width:768px){.layout{flex-direction:column}.nav{position:static;width:auto}}
  </style>
</head>
<body>
<div class="container">
  ${LOGO_HTML_TOP}

  <div class="header">
    <h1>~{workflow_name}</h1>
    <div class="meta">Version ~{version} • Run: ${RUN_DATE}</div>
  </div>

  <div class="layout">
    <nav class="nav">
      <strong>Sections</strong>
      <a href="#summary">Summary</a>
      <a href="#qc">Quality Control</a>
      <a href="#assembly">Assembly</a>
      <a href="#annotation">Annotation</a>
      <a href="#pangenome">Pangenome</a>
      <a href="#mlst">MLST</a>
      <a href="#variants">Variants</a>
      <a href="#amr">AMR</a>
      <a href="#mge">MGE</a>
      <a href="#virulence">Virulence</a>
      <a href="#blast">BLAST</a>
      <a href="#phylogeny">Phylogeny</a>
    </nav>

    <div class="content">
      <div id="summary" class="section">
        <h2>Analysis Summary</h2>
        <div class="report-card">
          <p>This report aggregates results from rMAP modules (QC, assembly/annotation, pangenome/phylogeny, AMR/MGE/virulence/BLAST, and phylogeny).</p>
        </div>
      </div>

      <div id="qc" class="section">
        <h2>Quality Control</h2>
EOF
    if (( ~{length(quality_reports)} )); then
      echo '        <div class="report-card"><h3>Per-sample QC reports</h3><ul>' >> "final_report/${REPORT_FILENAME}"
      idx=1
      for qf in ~{sep=' ' quality_reports}; do
        base="$(basename "$qf")"
        if [[ "${base,,}" == *combined* ]]; then continue; fi
        sid="$(extract_sample_id "$qf")"
        ext="${qf##*.}"
        unique="${sid}_qc.${ext}"
        cp -v "$qf" "final_report/assets/sections/$unique" || true
        echo "          <li><a href=\"assets/sections/$unique\" target=\"_blank\">QC Report ${idx} — ${sid}</a></li>" >> "final_report/${REPORT_FILENAME}"
        idx=$((idx+1))
      done
      echo '        </ul></div>' >> "final_report/${REPORT_FILENAME}"
    fi
    embed_html "Read trimming" "~{trimming_report_html}" "final_report/${REPORT_FILENAME}"
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # Assembly
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="assembly" class="section"><h2>Assembly</h2>
EOF
    embed_html "Assembly statistics" "~{assembly_stats_html}" "final_report/${REPORT_FILENAME}"
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # Annotation
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="annotation" class="section"><h2>Annotation</h2>
EOF
    embed_html "Annotation summary" "~{annotation_summary_html}" "final_report/${REPORT_FILENAME}"
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # Pangenome
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="pangenome" class="section"><h2>Pangenome</h2>
EOF
    if [[ -f "~{pangenome_report}" ]]; then
      embed_html "Pangenome report" "~{pangenome_report}" "final_report/${REPORT_FILENAME}"
    else
      embed_html "Pangenome report" "~{pangenome_report_html}" "final_report/${REPORT_FILENAME}"
    fi
    if [[ -n "${GH_BASENAME:-}" ]]; then
      cat >> "final_report/${REPORT_FILENAME}" <<EOF
        <div class="report-card">
          <h3>Gene Presence/Absence Heatmap</h3>
          <div class="image-center"><img src="assets/images/${GH_BASENAME}" alt="Gene heatmap" onerror="this.style.display='none'"></div>
        </div>
EOF
    fi
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # MLST
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="mlst" class="section"><h2>MLST</h2>
EOF
    embed_html "Combined MLST" "~{mlst_combined_html}" "final_report/${REPORT_FILENAME}"
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # Variants
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="variants" class="section"><h2>Variants</h2>
EOF
    embed_html "Variant summary" "~{variant_summary_html}" "final_report/${REPORT_FILENAME}"
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # AMR — per-sample links (skip "combined")
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="amr" class="section"><h2>Antimicrobial Resistance</h2>
EOF
    if (( ~{length(amr_reports)} )); then
      echo '        <div class="report-card"><h3>Per-sample AMR reports</h3><ul>' >> "final_report/${REPORT_FILENAME}"
      aidx=1
      for ar in ~{sep=' ' amr_reports}; do
        base="$(basename "$ar")"
        if [[ "${base,,}" == *combined* ]]; then continue; fi
        sid="$(extract_sample_id "$ar")"
        ext="${ar##*.}"
        unique="${sid}_amr.${ext}"
        cp -v "$ar" "final_report/assets/sections/$unique" || true
        echo "          <li><a href=\"assets/sections/$unique\" target=\"_blank\">AMR Report ${aidx} — ${sid}</a></li>" >> "final_report/${REPORT_FILENAME}"
        aidx=$((aidx+1))
      done
      echo '        </ul></div>' >> "final_report/${REPORT_FILENAME}"
    fi
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # MGE — per-sample links (skip "combined")
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="mge" class="section"><h2>Mobile Genetic Elements</h2>
EOF
    if (( ~{length(mge_reports)} )); then
      echo '        <div class="report-card"><h3>Per-sample MGE reports</h3><ul>' >> "final_report/${REPORT_FILENAME}"
      midx=1
      for mr in ~{sep=' ' mge_reports}; do
        base="$(basename "$mr")"
        if [[ "${base,,}" == *combined* ]]; then continue; fi
        sid="$(extract_sample_id "$mr")"
        ext="${mr##*.}"
        unique="${sid}_mge.${ext}"
        cp -v "$mr" "final_report/assets/sections/$unique" || true
        echo "          <li><a href=\"assets/sections/$unique\" target=\"_blank\">MGE Report ${midx} — ${sid}</a></li>" >> "final_report/${REPORT_FILENAME}"
        midx=$((midx+1))
      done
      echo '        </ul></div>' >> "final_report/${REPORT_FILENAME}"
    fi
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # VIRULENCE — per-sample links (skip "combined")
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="virulence" class="section"><h2>Virulence</h2>
EOF
    if (( ~{length(virulence_reports)} )); then
      echo '        <div class="report-card"><h3>Per-sample Virulence reports</h3><ul>' >> "final_report/${REPORT_FILENAME}"
      vidx=1
      for vf in ~{sep=' ' virulence_reports}; do
        base="$(basename "$vf")"
        if [[ "${base,,}" == *combined* ]]; then continue; fi
        sid="$(extract_sample_id "$vf")"
        ext="${vf##*.}"
        unique="${sid}_virulence.${ext}"
        cp -v "$vf" "final_report/assets/sections/$unique" || true
        echo "          <li><a href=\"assets/sections/$unique\" target=\"_blank\">Virulence Report ${vidx} — ${sid}</a></li>" >> "final_report/${REPORT_FILENAME}"
        vidx=$((vidx+1))
      done
      echo '        </ul></div>' >> "final_report/${REPORT_FILENAME}"
    fi
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # BLAST — per-sample links (skip "combined")
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="blast" class="section"><h2>BLAST</h2>
EOF
    if (( ~{length(blast_reports)} )); then
      echo '        <div class="report-card"><h3>Per-sample BLAST reports</h3><ul>' >> "final_report/${REPORT_FILENAME}"
      bidx=1
      for br in ~{sep=' ' blast_reports}; do
        base="$(basename "$br")"
        if [[ "${base,,}" == *combined* ]]; then continue; fi
        sid="$(extract_sample_id "$br")"
        ext="${br##*.}"
        unique="${sid}_blast.${ext}"
        cp -v "$br" "final_report/assets/sections/$unique" || true
        echo "          <li><a href=\"assets/sections/$unique\" target=\"_blank\">BLAST Report ${bidx} — ${sid}</a></li>" >> "final_report/${REPORT_FILENAME}"
        bidx=$((bidx+1))
      done
      echo '        </ul></div>' >> "final_report/${REPORT_FILENAME}"
    fi
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # Phylogeny
    cat >> "final_report/${REPORT_FILENAME}" <<'EOF'
      <div id="phylogeny" class="section"><h2>Phylogeny</h2>
EOF
    for img in ~{sep=' ' tree_images}; do
      if [[ -f "${img}" ]]; then
        b=$(basename "${img}")
        cat >> "final_report/${REPORT_FILENAME}" <<EOF
        <div class="report-card">
          <h3>${b}</h3>
          <div class="image-center"><img src="assets/images/${b}" alt="Phylogenetic tree" onerror="this.style.display='none'"></div>
        </div>
EOF
      fi
    done
    echo "      </div>" >> "final_report/${REPORT_FILENAME}"

    # Footer (close HTML)
    cat >> "final_report/${REPORT_FILENAME}" <<EOF
    </div> <!-- content -->
  </div> <!-- layout -->

  ${LOGO_HTML_BOTTOM}
  <div class="footer-banner">FOOTER_SENTENCE_PLACEHOLDER</div>
</div>
</body></html>
EOF

    # Extra safety: if any literal ${RUN_DATE} survived (e.g., someone re-quoted a heredoc),
    perl -0777 -pe 's/\$\{RUN_DATE\}/$ENV{RUN_DATE}/g' -i "final_report/${REPORT_FILENAME}" || true

    FS="~{footer_sentence}"
    if [[ -z "${FS// }" ]]; then
      FS="Analysis performed using rMAP version ~{version}"
    else
      FS="${FS//\{RUN_DATE\}/${RUN_DATE}}"
      FS="${FS//\{VERSION\}/~{version}}"
    fi
    sed -i "s|FOOTER_SENTENCE_PLACEHOLDER|${FS}|" "final_report/${REPORT_FILENAME}"

    echo "Creating final_report.tgz..."
    tar -czvf final_report.tgz final_report || cp -v default_report.tgz final_report.tgz

    [[ -f "final_report.tgz" && -f "final_report/${REPORT_FILENAME}" ]] || { echo "Error: outputs missing" >&2; exit 1; }
    echo "=== Report generation complete $(date) ==="
    ls -lh final_report.tgz final_report/${REPORT_FILENAME} || true
  >>>

  runtime {
    docker: "ubuntu:20.04"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    dockerVolumeMode: "delegated"
    continueOnReturnCode: true
    maxRetries: 3
    preemptible: 0
    bootDiskSizeGb: 20
    bootDiskType: "pd-ssd"
  }

  output {
    File final_report_html = "final_report/rMAP_~{version}_final_report.html"
    File final_report_tgz  = "final_report.tgz"
    File skip_log          = "skip_report.log"
    File? generation_log   = "report_generation.log"
  }

  parameter_meta {
    trimming_report_html: "HTML report from read trimming"
    assembly_stats_html:  "HTML with assembly statistics"
    annotation_summary_html: "HTML with annotation summary"
    pangenome_report_html: "HTML pangenome report"
    gene_heatmap_png:     "PNG image of gene presence/absence heatmap"
    mlst_combined_html:   "Combined MLST results in HTML"
    variant_summary_html: "HTML variant summary"

    amr_reports:          "Array of per-sample AMR HTML reports"
    mge_reports:          "Array of per-sample MGE HTML reports"
    virulence_reports:    "Array of per-sample virulence HTML reports"
    quality_reports:      "Array of per-sample QC HTML reports"
    blast_reports:        "Array of per-sample BLAST HTML reports"

    pangenome_report:     "Direct pangenome report file"
    gene_heatmap:         "Direct gene heatmap image file"
    tree_images:          "Array of phylogenetic tree images"
    workflow_name:        "Name of the workflow for report title"
    version:              "Version number to display in report"
    footer_sentence:      "Custom footer text (supports {RUN_DATE} and {VERSION} placeholders)"
    skip:                 "Whether to skip report generation"
    logo_file:            "Optional logo file to use instead of downloading"
  }
}
