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
    Int max_cpus = 8
    Int max_memory_gb = 12
    Int min_assembly_quality = 50
    Int min_read_length = 50
    Int min_mapping_quality = 20
    Int tree_image_width = 1200
    String tree_image_format = "png"
    Int tree_font_size = 8

    # Derived values with adjusted resources
    Int cpu_4 = if (max_cpus < 4) then max_cpus else 4
    Int cpu_8 = if (max_cpus < 8) then max_cpus else max_cpus
    Int cpu_16 = max_cpus
    Int cpu_2 = if (max_cpus < 2) then max_cpus else 2
    Int mem_8 = if (max_memory_gb < 8) then max_memory_gb else 8
    Int mem_12 = max_memory_gb
  }

  meta {
    workflow_timeout: "168 hours"
    workflow_heartbeat_interval: "10 minutes"
    workflow_heartbeat_ttl: "30 minutes"
    allowNestedInputs: true
    maxRetries: 3
    continueOnReturnCode: [0, 1]
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
      memory_gb = mem_12,
      min_quality = min_assembly_quality
  }

  call ANNOTATION {
    input:
      assembly_output = ASSEMBLY.assembly_output,
      do_annotation = do_annotation,
      cpu = cpu_8
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
      min_quality = min_mapping_quality
  }

  call PANGENOME {
    input:
      annotation_input = ANNOTATION.annotation_output,
      do_pangenome = do_pangenome,
      cpu = cpu_8,
      memory_gb = mem_12,
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
      contig_fastas = ASSEMBLY.assembly_output,  # Direct array
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
        memory_gb = mem_12,
        tree_prefix = "core_genes",
        bootstrap_replicates = 100
    }

    call ACCESSORY_PHYLOGENY {
      input:
        alignment = PANGENOME.accessory_binary,
        do_phylogeny = true,
        model = phylogeny_model,
        cpu = cpu_8,
        tree_prefix = "accessory_genes",
        bootstrap_replicates = 100
    }

    call TREE_VISUALIZATION as CORE_TREE {
      input:
        input_tree  = CORE_PHYLOGENY.phylogeny_tree,
        width       = tree_image_width,
        image_format= tree_image_format,
        font_size   = tree_font_size,
        tree_title  = "Core Genes Phylogenetic Tree"
    }

    call TREE_VISUALIZATION as ACCESSORY_TREE {
      input:
        input_tree  = ACCESSORY_PHYLOGENY.phylogeny_tree,
        width       = tree_image_width,
        image_format= tree_image_format,
        font_size   = tree_font_size,
        tree_title  = "Accessory Genes Phylogenetic Tree"
    }
  }

  # --- locals ---
  Array[File] blast_reports_local = BLAST_ANALYSIS.blast_reports
  Array[File] tree_images_local = if do_phylogeny then select_all([CORE_TREE.final_image, ACCESSORY_TREE.final_image]) else []
  File?       pangenome_report_html_local = PANGENOME.pangenome_report
  File?       gene_heatmap_png_local      = PANGENOME.gene_heatmap

  if (do_reporting) {
    call MERGE_REPORTS {
      input:
        trimming_report_html   = TRIMMING.trimming_report,
        fastqc_summary_html    = QUALITY_CONTROL.fastqc_summary_html,
        assembly_stats_html    = ASSEMBLY.assembly_stats,
        annotation_summary_html= ANNOTATION.annotation_summary,
        pangenome_report_html  = pangenome_report_html_local,
        gene_heatmap_png       = gene_heatmap_png_local,
        pangenome_report       = PANGENOME.pangenome_report,
        gene_heatmap           = PANGENOME.gene_heatmap,
        mlst_combined_html     = MLST.combined_html,
        variant_summary_html   = VARIANT_CALLING.variant_summary,
        amr_combined_html      = AMR_PROFILING.combined_html,
        mge_combined_html      = MGE_ANALYSIS.html_report,
        virulence_combined_html= VIRULENCE_ANALYSIS.combined_html,
        blast_reports          = blast_reports_local,
        tree_images            = tree_images_local,
        workflow_name          = "rMAP Analysis Pipeline",
        version                = "2.0",
        footer_sentence        = "This report was generated by rMAP 2.0 pipeline"
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

    File core_phylogeny_output      = select_first([CORE_PHYLOGENY.phylogeny_tree, "default_core_tree.nwk"])
    File accessory_phylogeny_output = select_first([ACCESSORY_PHYLOGENY.phylogeny_tree, "default_accessory_tree.nwk"])

    # New: explicit rendered images + logs (one per tree)
    File? core_tree_image           = CORE_TREE.final_image
    File? accessory_tree_image      = ACCESSORY_TREE.final_image
    File? core_tree_render_log      = CORE_TREE.render_log
    File? accessory_tree_render_log = ACCESSORY_TREE.render_log

    Array[File] amr_output   = AMR_PROFILING.amr_outputs
    File        combined_amr = AMR_PROFILING.combined_amr

    File        plasmid_report  = MGE_ANALYSIS.plasmid_report
    Array[File] plasmid_results = MGE_ANALYSIS.sample_reports

    Array[File] blast_results = BLAST_ANALYSIS.blast_results
    Array[File] blast_top10   = BLAST_ANALYSIS.blast_top10
    Array[File] blast_logs    = BLAST_ANALYSIS.blast_logs

    Array[File] virulence_reports         = VIRULENCE_ANALYSIS.virulence_reports
    File        combined_virulence_report = VIRULENCE_ANALYSIS.combined_report

    # Final report
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

    if [ "~{skip}" == "true" ]; then
      echo "Skipping trimming process as requested" > trimming_skipped.log
      echo "Creating empty output files for portability" >> trimming_skipped.log

      mkdir -p trimmed

      # Create empty output files for each input pair
      counter=0
      R1=""
      for file in ~{sep=' ' input_reads}; do
        if [ $((counter % 2)) -eq 0 ]; then
          R1="$file"
        else
          R2="$file"
          R1_base=$(basename "$R1")
          sample_name=$(echo "$R1_base" | sed -e 's/[._][Rr]1[._].*//' -e 's/[._]1[._].*//')

          # Create empty output files
          touch "trimmed/${sample_name}_1.trim.fastq.gz"
          touch "trimmed/${sample_name}_2.trim.fastq.gz"
          echo "Created empty output files for sample: $sample_name" >> trimming_skipped.log
        fi
        counter=$((counter + 1))
      done

      # Create empty statistics files
      echo "Sample,Input_Pairs,Both_Surviving,Forward_Only,Reverse_Only,Dropped" > trimming_stats.csv
      echo "0,0,0,0,0" >> trimming_stats.csv

      # Create minimal HTML report
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

      # Initialize statistics file
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

          # Create temporary log file for this sample
          sample_log="trimmed/${sample_name}.log"
          touch "$sample_log"

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
              MINLEN:~{min_length} > "$sample_log" 2>&1 || {
                echo "WARNING: Trimmomatic failed for sample $sample_name" >> trimming.log
                # Create empty output files to allow pipeline to continue
                touch "trimmed/${sample_name}_1.trim.fastq.gz" "trimmed/${sample_name}_2.trim.fastq.gz"
              }
          )

          # Extract statistics from the sample log
          stats_line=$(grep "Input Read Pairs" "$sample_log" || true)
          if [[ -n "${stats_line:-}" ]]; then
            # Extract numbers using awk
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

          # Append sample log to main log
          cat "$sample_log" >> trimming.log
          rm "$sample_log"

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

      # Generate HTML report
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
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Input Pairs</th>
                <th>Both Surviving</th>
                <th>Forward Only</th>
                <th>Reverse Only</th>
                <th>Dropped</th>
                <th>Survival Rate</th>
            </tr>
        </thead>
        <tbody>
EOF

      # Add data rows to HTML
      while IFS=',' read -r sample input both forward reverse dropped; do
        # Skip header line
        if [ "$sample" == "Sample" ]; then
          continue
        fi

        # Calculate survival rate
        if [ "$input" -gt 0 ]; then
          total_surviving=$((both + forward + reverse))
          survival_rate=$(awk -v ts="$total_surviving" -v i="$input" 'BEGIN {printf "%.2f", ts*100/i}')
        else
          survival_rate="0.00"
        fi

        # Determine CSS class
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
        echo "                <td>$dropped</td>" >> trimming_report.html
        echo "                <td class=\"$class\">${survival_rate}%</td>" >> trimming_report.html
        echo "            </tr>" >> trimming_report.html
      done < trimming_stats.csv

      # Close HTML
      cat >> trimming_report.html <<EOF
        </tbody>
    </table>
</body>
</html>
EOF

      echo "Trimming completed successfully at $(date)" >> trimming.log
      echo "Output files created:" >> trimming.log
      ls -lh trimmed/* >> trimming.log

      # Ensure skip marker exists even when not skipping, so outputs never break
      echo "Trimming executed (not skipped) on $(date)" > trimming_skipped.log
    else
      echo "Trimming skipped by user request" > trimming_skipped.txt
      # Also ensure skip marker exists for outputs
      echo "Trimming skipped on $(date)" > trimming_skipped.log
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
    File? trimming_stats = if do_trimming then "trimming_stats.csv" else "trimming_skipped.txt"
    File? trimming_report = if do_trimming then "trimming_report.html" else "trimming_skipped.txt"
    File skip_log = "trimming_skipped.log"
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
    <title>FastQC Summary Report</title>
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
    <h1>FastQC Summary Report</h1>
    <p>Generated at $(date)</p>
    <table>
        <thead>
            <tr>
                <th>Sample</th>
                <th>Report</th>
                <th>Data</th>
            </tr>
        </thead>
        <tbody>
EOF

          # Add entries for each FastQC report
          for f in qc_reports/*_fastqc.html; do
            sample=$(basename "$f" | sed 's/_fastqc\.html//')
            echo "            <tr>" >> qc_reports/fastqc_summary.html
            echo "                <td>$sample</td>" >> qc_reports/fastqc_summary.html
            echo "                <td><a href=\"$(basename "$f")\">HTML Report</a></td>" >> qc_reports/fastqc_summary.html
            zip_file="${f%.html}.zip"
            echo "                <td><a href=\"$(basename "$zip_file")\">Raw Data</a></td>" >> qc_reports/fastqc_summary.html
            echo "            </tr>" >> qc_reports/fastqc_summary.html
          done

          cat >> qc_reports/fastqc_summary.html <<EOF
        </tbody>
    </table>
    <p>Note: Full MultiQC report not available. Install MultiQC for more comprehensive analysis.</p>
</body>
</html>
EOF
        else
          # Run MultiQC if available
          multiqc qc_reports -o qc_reports --force 2>> qc.log || {
            echo "ERROR: MultiQC failed" >> qc.log
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
    Array[File] quality_reports = if do_quality_control then glob("qc_reports/*.{html,zip,txt}") else glob("qc_reports/skipped.txt")
    File? qc_log = if do_quality_control then "qc.log" else "qc_reports/skipped.txt"
    File? fastqc_summary_html = if do_quality_control then "qc_reports/summary.html" else "qc_reports/skipped.txt"
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
    Int memory_gb = 8
    Int min_quality = 50
    Boolean run_quast = false  # Changed to false since we can't guarantee QUAST availability
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
        sample_name=$(basename "$R1" | sed 's/_1.fastq.gz//; s/_1.trim.fastq.gz//')
        outdir="~{output_dir}/megahit_${sample_name}"
        contig_file="~{output_dir}/${sample_name}_contigs.fa"

        echo "Assembling $sample_name (R1: $R1, R2: $R2)" >> assembly.log

        # Run assembly in a subshell
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
        sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
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

    # Make a blank heatmap image if we can (optional)
    if command -v convert >/dev/null 2>&1; then
      convert -size ~{heatmap_width}x~{heatmap_height} xc:white \
              -gravity center -pointsize 24 -annotate 0 "No Data" \
              final_output/gene_presence_heatmap.png || true
    else
      : > final_output/gene_presence_heatmap.png || true
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

        # Try to draw a real heatmap if R is available; otherwise keep placeholder
        if command -v Rscript >/dev/null 2>&1; then
          cat > final_output/generate_visualizations.R <<'RS'
suppressWarnings(suppressMessages({
  try({
    library(pheatmap); library(RColorBrewer)
    df <- read.csv("gene_presence_absence.csv", header=TRUE, check.names=FALSE)
    # Heuristic for sample columns
    sample_cols <- setdiff(colnames(df), colnames(df)[1:14])
    if (length(sample_cols) == 0) sample_cols <- colnames(df)
    mat <- as.matrix(df[, sample_cols, drop=FALSE])
    mat[mat != ""] <- 1; mat[mat == ""] <- 0
    mode(mat) <- "numeric"
    png("gene_presence_heatmap.png", width=~{heatmap_width}, height=~{heatmap_height})
    pheatmap(mat, color = colorRampPalette(brewer.pal(9,"Blues"))(2),
             cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
             main="Gene Presence/Absence")
    dev.off()
  }, silent=TRUE)
})
RS
          ( cd final_output && Rscript generate_visualizations.R ) || true
        fi

        # Stamp success into the HTML
        sed -i 's/Initialized./Roary completed successfully./' final_output/pangenome_report.html || true
      else
        sed -i 's/Initialized./Roary was invoked but did not complete successfully; placeholders shown./' final_output/pangenome_report.html || true
      fi
    else
      sed -i 's/Initialized./Roary not available in container; placeholders shown./' final_output/pangenome_report.html || true
    fi
  >>>

  runtime {
    docker: "staphb/roary:3.13.0"
    memory: "~{memory_gb}G"
    cpu: cpu
    disks: "local-disk 100 HDD"
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
    File phylogeny_tree = "phylogeny_results/core_~{tree_prefix}.nwk"
    # expose both logs as optional; one or both may exist
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
    set -euo pipefail

    # Initialize directories and logging
    mkdir -p phylogeny_results
    echo "Accessory Phylogeny Analysis Log - $(date)" > phylogeny.log
    echo "=========================================" >> phylogeny.log
    echo "Runtime Parameters:" >> phylogeny.log
    echo "- do_phylogeny: ~{do_phylogeny}" >> phylogeny.log
    echo "- tree_prefix: ~{tree_prefix}" >> phylogeny.log
    echo "- model: ~{model}" >> phylogeny.log
    echo "- cpu: ~{cpu}" >> phylogeny.log
    echo "- bootstrap_replicates: ~{bootstrap_replicates}" >> phylogeny.log
    echo "- memory_gb: ~{memory_gb}" >> phylogeny.log
    echo "=========================================" >> phylogeny.log

    # Create skipped.txt file for cases where analysis is skipped
    echo "Analysis skipped" > skipped.txt

    # Ensure error log file exists early so outputs can bind cleanly
    : > "phylogeny_results/accessory_error.log"

    # Skip condition 1: User explicitly disabled phylogeny
    if [ "~{do_phylogeny}" != "true" ]; then
      echo "Accessory phylogeny disabled by user parameter" >> phylogeny.log
      echo "(ACCESSORY_PHYLOGENY_DISABLED);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Analysis skipped by user request" > "phylogeny_results/accessory_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      # Create a canonical log name even in skip path
      cp "phylogeny_results/accessory_~{tree_prefix}.log" "phylogeny_results/accessory_~{tree_prefix}.log"
      exit 0
    fi

    # Skip condition 2: Alignment file not provided
    if [ ! -f "~{alignment}" ]; then
      echo "ERROR: Accessory alignment file not found at path: ~{alignment}" >> phylogeny.log
      echo "(MISSING_ACCESSORY_ALIGNMENT);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Accessory alignment file missing" > "phylogeny_results/accessory_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      exit 0
    fi

    # Skip condition 3: Alignment file exists but is empty
    if [ ! -s "~{alignment}" ]; then
      echo "ERROR: Accessory alignment file is empty: ~{alignment}" >> phylogeny.log
      echo "(EMPTY_ACCESSORY_ALIGNMENT);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Accessory alignment file empty" > "phylogeny_results/accessory_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      exit 0
    fi

    # Validate alignment content
    echo "Validating accessory alignment file..." >> phylogeny.log
    seq_count=$(grep -c '^>' "~{alignment}" || echo 0)
    echo "Found $seq_count sequences in accessory alignment" >> phylogeny.log

    # Skip condition 4: Insufficient sequences
    if [ "$seq_count" -lt 4 ]; then
      echo "ERROR: Insufficient sequences ($seq_count) for accessory phylogeny (minimum 4 required)" >> phylogeny.log
      echo "(INSUFFICIENT_ACCESSORY_SEQUENCES_$seq_count);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      echo "Only $seq_count accessory sequences found" > "phylogeny_results/accessory_~{tree_prefix}.log"
      echo "System Info:" >> phylogeny.log
      free -h >> phylogeny.log
      exit 0
    fi

    # Run phylogenetic analysis
    echo "Starting accessory phylogenetic analysis with FastTree..." >> phylogeny.log
    echo "System memory information:" >> phylogeny.log
    free -h >> phylogeny.log

    set +e
    ulimit -v $((~{memory_gb} * 1024 * 1024))

    echo "Command: FastTree ~{model} -gamma -quiet -boot ~{bootstrap_replicates} \\" >> phylogeny.log
    echo "  -log phylogeny_results/accessory_~{tree_prefix}.log \\" >> phylogeny.log
    echo "  < ~{alignment} > phylogeny_results/accessory_~{tree_prefix}.nwk" >> phylogeny.log

    FastTree ~{model} \
      -gamma \
      -quiet \
      -boot ~{bootstrap_replicates} \
      -log "phylogeny_results/accessory_~{tree_prefix}.log" \
      < "~{alignment}" > "phylogeny_results/accessory_~{tree_prefix}.nwk" 2> "phylogeny_results/accessory_error.log"
    exit_code=$?
    set -e

    # Handle FastTree failure
    if [ $exit_code -ne 0 ]; then
      echo "WARNING: FastTree exited with code $exit_code" >> phylogeny.log
      cat "phylogeny_results/accessory_error.log" >> phylogeny.log

      # Attempt reduced bootstrap replicates if memory issue
      if grep -qi "oom" phylogeny.log || grep -qi "killed" phylogeny.log; then
        echo "Attempting with reduced bootstrap replicates (50)..." >> phylogeny.log
        FastTree ~{model} \
          -gamma \
          -quiet \
          -boot 50 \
          -log "phylogeny_results/accessory_~{tree_prefix}_reduced.log" \
          < "~{alignment}" > "phylogeny_results/accessory_~{tree_prefix}.nwk" 2>> phylogeny.log || {
            echo "Fallback failed, generating minimal tree" >> phylogeny.log
            echo "(ACCESSORY_TREE_FAILED);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
          }
        # Normalize log name if reduced run created a different log file
        if [ -f "phylogeny_results/accessory_~{tree_prefix}_reduced.log" ]; then
          cp "phylogeny_results/accessory_~{tree_prefix}_reduced.log" "phylogeny_results/accessory_~{tree_prefix}.log"
        fi
      else
        echo "Generating minimal tree after failure" >> phylogeny.log
        echo "(ACCESSORY_TREE_FAILED);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
      fi
    fi

    # Final validation of output
    if [ ! -s "phylogeny_results/accessory_~{tree_prefix}.nwk" ]; then
      echo "ERROR: Accessory tree file is empty, generating minimal tree" >> phylogeny.log
      echo "(EMPTY_ACCESSORY_OUTPUT);" > "phylogeny_results/accessory_~{tree_prefix}.nwk"
    fi

    # Guarantee canonical log file exists (even if none were produced earlier)
    if [ ! -s "phylogeny_results/accessory_~{tree_prefix}.log" ]; then
      : > "phylogeny_results/accessory_~{tree_prefix}.log"
    fi

    echo "Accessory phylogenetic analysis completed at $(date)" >> phylogeny.log
    echo "Final output files:" >> phylogeny.log
    ls -lh phylogeny_results/* >> phylogeny.log
  >>>

  runtime {
    docker: "staphb/fasttree:2.1.11"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
    timeout: "12 hours"
  }

  output {
    File phylogeny_tree = "phylogeny_results/accessory_~{tree_prefix}.nwk"
    File? phylogeny_log_reduced = "phylogeny_results/accessory_~{tree_prefix}_reduced.log"
    File? phylogeny_log = "phylogeny_results/accessory_~{tree_prefix}.log"
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
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
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
    echo "=====================================" >> variant.log

    # Function to ensure HTML output exists
    ensure_output() {
      local status=$1
      mkdir -p variants

      if [ "$status" == "success" ]; then
        # Create minimal HTML if final version wasn't generated
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
        # Create error HTML
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

    # Process samples
    RESULTS_TEMP=$(mktemp) || { echo "ERROR: Failed to create temp file" >> variant.log; exit 1; }
    echo -e "Sample\tStatus\tSNPs\tINDELs\tTotal" > "$RESULTS_TEMP"
    processed_samples=0

    for ((i=0; i<${#files[@]}; i+=2)); do
      R1="${files[i]}"
      R2="${files[i+1]}"
      SAMPLE_NAME=$(basename "$R1" | sed 's/[._][Rr]1.*//; s/[._]1.*//; s/[._][12].*//')
      SAMPLE_DIR="variants/$SAMPLE_NAME"

      mkdir -p "$SAMPLE_DIR" || {
        echo "ERROR: Failed to create directory for $SAMPLE_NAME" >> variant.log
        echo -e "$SAMPLE_NAME\tfailed\t0\t0\t0" >> "$RESULTS_TEMP"
        continue
      }

      echo "Processing $SAMPLE_NAME" >> variant.log
      SNP_COUNT=0
      INDEL_COUNT=0
      TOTAL=0
      STATUS="success"

      # Handle reference type
      if [ "~{reference_type}" == "genbank" ]; then
        ref_path="$SAMPLE_DIR/reference.gbk"
      else
        ref_path="$SAMPLE_DIR/reference.fasta"
      fi

      if ! cp "~{reference_genome}" "$ref_path"; then
        echo "ERROR: Failed to copy reference for $SAMPLE_NAME" >> variant.log
        STATUS="failed"
      else
        # Run variant calling
        set +e
        snippy --cpus ~{cpu} --minqual ~{min_quality} \
               --ref "$ref_path" --R1 "$R1" --R2 "$R2" \
               --outdir "$SAMPLE_DIR" --prefix "$SAMPLE_NAME" --force 2>> variant.log
        snippy_exit=$?
        set -e

        if [ $snippy_exit -ne 0 ] || [ ! -s "$SAMPLE_DIR/$SAMPLE_NAME.vcf" ]; then
          STATUS="failed"
          echo "WARNING: snippy failed for $SAMPLE_NAME (exit $snippy_exit)" >> variant.log
        fi
      fi

      # Count variants if successful
      if [ "$STATUS" == "success" ]; then
        SNP_COUNT=$(grep -c 'TYPE=snp' "$SAMPLE_DIR/$SAMPLE_NAME.vcf" 2>/dev/null || echo 0)
        INDEL_COUNT=$(awk '/TYPE=indel/ {count++} END {print count+0}' "$SAMPLE_DIR/$SAMPLE_NAME.vcf" 2>/dev/null)
        TOTAL=$((SNP_COUNT + INDEL_COUNT))
        processed_samples=$((processed_samples + 1))
      fi

      echo -e "$SAMPLE_NAME\t$STATUS\t$SNP_COUNT\t$INDEL_COUNT\t$TOTAL" >> "$RESULTS_TEMP"
      echo "Processed $SAMPLE_NAME: Status=$STATUS, SNPs=$SNP_COUNT, INDELs=$INDEL_COUNT, Total=$TOTAL" >> variant.log
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
      </tr>
    </thead>
    <tbody>
EOF

    while IFS=$'\t' read -r sample status snps indels total; do
      if [ "$sample" == "Sample" ]; then continue; fi
      if [ "$status" == "failed" ]; then
        echo "<tr class=\"error\"><td>$sample</td><td>Failed</td><td>-</td><td>-</td><td>-</td></tr>" >> variants/variant_summary.html
      else
        if [ $total -eq 0 ]; then
          row_class="success"
        elif [ $total -lt 100 ]; then
          row_class="warning"
        else
          row_class="error"
        fi
        echo "<tr><td>$sample</td><td>Success</td><td>$snps</td><td>$indels</td><td class=\"$row_class\">$total</td></tr>" >> variants/variant_summary.html
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
  >>>

  runtime {
    docker: "quay.io/biocontainers/snippy:4.6.0--hdfd78af_1"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 100 HDD"
    continueOnReturnCode: true
    preemptible: 2
    timeout: "12 hours"
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
  }

  command <<<
    set -euo pipefail
    set -x

    mkdir -p amr_results html_results
    echo "AMR Profiling Analysis Log - $(date)" > amr.log
    echo "===================================" >> amr.log
    echo "Runtime Parameters:" >> amr.log
    echo "- do_amr_profiling: ~{do_amr_profiling}" >> amr.log
    echo "- use_local_db: ~{use_local_db}" >> amr.log
    echo "- minid: ~{minid}" >> amr.log
    echo "- mincov: ~{mincov}" >> amr.log
    echo "- cpu: ~{cpu}" >> amr.log
    echo "===================================" >> amr.log

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

    processed_samples=0
    for asm_file in ~{sep=' ' assembly_output}; do
      [ ! -f "$asm_file" ] && continue
      [ ! -s "$asm_file" ] && continue

      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//' | sed 's/_contigs//')
      output_file="amr_results/${sample_name}_amr.tsv"

      abricate \
        --db $db_to_use \
        --minid ~{minid} \
        --mincov ~{mincov} \
        --threads ~{cpu} \
        --csv \
        "$asm_file" > "${output_file}.tmp" 2>> amr.log || {
          echo "WARNING: AMR profiling failed for $sample_name" >> amr.log
          echo -e "SAMPLE\tCONTIG\tGENE\t%COVERAGE\t%IDENTITY\tPRODUCT\tRESISTANCE" > "$output_file"
          echo -e "$sample_name\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR" >> "$output_file"
          continue
        }

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

      processed_samples=$((processed_samples + 1))
      rm -f "${output_file}.tmp"
    done
  >>>

  runtime {
    docker: "staphb/abricate:1.0.0"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
    timeout: "6 hours"
  }

  output {
    Array[File] amr_outputs = if do_amr_profiling then glob("amr_results/*_amr.tsv") else ["amr_results/skipped.txt"]
    File combined_amr = if do_amr_profiling then "amr_results/combined_amr.tsv" else "amr_results/skipped.txt"
    File combined_html = if do_amr_profiling then "html_results/combined_amr.html" else "html_results/skipped.html"
    File amr_log = "amr.log"
  }
}

task MGE_ANALYSIS {
  input {
    Array[File]? assembly_output
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

    # Skip condition 1: User explicitly disabled MGE analysis
    if [ "~{do_mge_analysis}" != "true" ]; then
      echo "MGE analysis disabled by user parameter" >> mge.log
      echo "MGE analysis skipped by user request" > mge_results/skipped.txt
      echo "<h1>MGE analysis skipped by user request</h1>" > html_results/skipped.html
      exit 0
    fi

    # Skip condition 2: No input files provided
    if [ -z "~{sep=' ' assembly_output}" ]; then
      echo "ERROR: No assembly files provided for MGE analysis" >> mge.log
      echo "NO_INPUT_FILES" > mge_results/skipped.txt
      echo "<h1>MGE analysis skipped - no input files provided</h1>" > html_results/skipped.html
      exit 0
    fi

    # Verify input files
    echo "Input files verification:" >> mge.log
    valid_files=0
    for f in ~{sep=' ' assembly_output}; do
      if [ ! -f "$f" ]; then
        echo "WARNING: Input file not found: $f" >> mge.log
      else
        size=$(wc -c < "$f")
        if [ $size -gt 0 ]; then
          echo "- Valid input file: $f ($size bytes)" >> mge.log
          valid_files=$((valid_files + 1))
        else
          echo "WARNING: Empty input file: $f" >> mge.log
        fi
      fi
    done

    # Skip condition 3: No valid input files
    if [ $valid_files -eq 0 ]; then
      echo "ERROR: No valid input files available" >> mge.log
      echo "NO_VALID_INPUTS" > mge_results/skipped.txt
      echo "<h1>MGE analysis skipped - no valid input files</h1>" > html_results/skipped.html
      exit 0
    fi

    # HTML template components
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
        document.getElementById("samples-processed").textContent = "'$valid_files'";

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
    else
      if ! abricate --list | grep -q "plasmidfinder"; then
        if ! abricate --setupdb --db plasmidfinder >> mge.log 2>&1; then
          echo "ERROR: Failed to setup plasmidfinder database" >> mge.log
          echo "DB_SETUP_FAILED" > mge_results/skipped.txt
          echo "<h1>MGE analysis skipped - database setup failed</h1>" > html_results/skipped.html
          exit 0
        fi
      fi
    fi

    # Process each assembly file
    processed_samples=0
    for asm_file in ~{sep=' ' assembly_output}; do
      [ ! -f "$asm_file" ] && continue
      [ ! -s "$asm_file" ] && continue

      sample_name=$(basename "$asm_file" | sed 's/\..*//')
      output_tsv="mge_results/${sample_name}_mge.tsv"
      output_html="html_results/${sample_name}_mge.html"

      echo "Processing $sample_name" >> mge.log

      # Run MGE analysis with error handling
      set +e
      abricate \
        --db $db_to_use \
        --mincov 80 \
        --minid 90 \
        --threads ~{cpu} \
        --nopath \
        "$asm_file" > "${output_tsv}.tmp" 2>> mge.log
      status=$?
      set -e

      if [ $status -ne 0 ] || [ ! -s "${output_tsv}.tmp" ]; then
        echo "WARNING: MGE analysis failed for $sample_name" >> mge.log
        echo -e "Sample\tMGE_Type\tGene\tProduct\t%Coverage\t%Identity\tAccession" > "$output_tsv"
        echo -e "$sample_name\tERROR\tERROR\tERROR\tERROR\tERROR\tERROR" >> "$output_tsv"
      else
        # Process successful output
        awk -v sample="$sample_name" '
        BEGIN {
          OFS="\t";
          print "Sample\tMGE_Type\tGene\tProduct\t%Coverage\t%Identity\tAccession";
        }
        NR==1 { next }
        {
          print sample, "Plasmid", $5, $13, $9, $10, $12;
        }' "${output_tsv}.tmp" > "$output_tsv"
        processed_samples=$((processed_samples + 1))
      fi
      rm -f "${output_tsv}.tmp"

      # Generate HTML report for this sample
      {
        echo "$HTML_HEADER"

        if grep -q "ERROR" "$output_tsv"; then
          echo "<tr class=\"error-row\">"
          echo "<td colspan=\"7\">MGE analysis failed for $sample_name</td>"
          echo "</tr>"
        else
          tail -n +2 "$output_tsv" | while IFS=$'\t' read -r _ mge_type gene product coverage identity accession; do
            echo "<tr>"
            echo "<td>$sample_name</td>"
            echo "<td>$mge_type</td>"
            echo "<td class=\"mge\">$gene</td>"
            echo "<td>$product</td>"
            echo "<td>$coverage</td>"
            echo "<td>$identity</td>"
            echo "<td>$accession</td>"
            echo "</tr>"
          done
        fi

        echo "$HTML_FOOTER"
      } > "$output_html"
    done

    # Generate combined reports if we processed any samples
    if [ $processed_samples -gt 0 ]; then
      echo -e "Sample\tMGE_Type\tGene\tProduct\t%Coverage\t%Identity\tAccession" > mge_results/combined_mge.tsv
      for tsv_file in mge_results/*_mge.tsv; do
        tail -n +2 "$tsv_file" >> mge_results/combined_mge.tsv
      done

      {
        echo "$HTML_HEADER"

        while IFS=$'\t' read -r sample mge_type gene product coverage identity accession; do
          if [ "$sample" == "Sample" ]; then continue; fi

          if [[ "$gene" == *"ERROR"* ]]; then
            echo "<tr class=\"error-row\">"
            echo "<td>$sample</td>"
            echo "<td colspan=\"6\">MGE analysis failed</td>"
            echo "</tr>"
          else
            echo "<tr>"
            echo "<td>$sample</td>"
            echo "<td>$mge_type</td>"
            echo "<td class=\"mge\">$gene</td>"
            echo "<td>$product</td>"
            echo "<td>$coverage</td>"
            echo "<td>$identity</td>"
            echo "<td>$accession</td>"
            echo "</tr>"
          fi
        done < mge_results/combined_mge.tsv

        echo "$HTML_FOOTER"
      } > html_results/combined_mge.html
    else
      echo "ERROR: No samples processed successfully" >> mge.log
      echo "NO_SUCCESSFUL_SAMPLES" > mge_results/combined_mge.tsv
      echo "<h1>No samples processed successfully</h1>" > html_results/combined_mge.html
    fi

    echo "MGE analysis completed at $(date)" >> mge.log
    echo "Processed $processed_samples samples successfully" >> mge.log
    echo "Output files:" >> mge.log
    ls -lh mge_results/* html_results/* >> mge.log 2>&1 || true
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk 75 HDD"
    continueOnReturnCode: true
    timeout: "12 hours"
  }

  output {
    Array[File] sample_reports = if do_mge_analysis then glob("mge_results/*_mge.tsv") else ["mge_results/skipped.txt"]
    File plasmid_report = if do_mge_analysis then "mge_results/combined_mge.tsv" else "mge_results/skipped.txt"
    File html_report = if do_mge_analysis then "html_results/combined_mge.html" else "html_results/skipped.html"
    File mge_log = "mge.log"
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
      echo "Creating empty output files for portability" >> virulence.log

      mkdir -p virulence_results
      mkdir -p html_results

      # Create empty output files
      touch virulence_results/combined_virulence.tsv
      touch html_results/combined_virulence.html

      # Create minimal TSV header
      echo -e "Sample\tVirulence_Factor\tProduct\t%Coverage\t%Identity\tRisk_Level" > virulence_results/combined_virulence.tsv

      # Create minimal HTML
      echo '<!DOCTYPE html>
<html>
<head>
    <title>Virulence Factor Analysis</title>
</head>
<body>
    <h1>Analysis Skipped</h1>
    <p>Virulence analysis was skipped as requested in the workflow parameters.</p>
</body>
</html>' > html_results/combined_virulence.html

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

    for asm_file in ~{sep=' ' assembly_output}; do
      sample_name=$(basename "$asm_file" | sed 's/\.[^.]*$//')
      output_tsv="virulence_results/${sample_name}_virulence.tsv"
      output_html="html_results/${sample_name}_virulence.html"

      echo "Processing $sample_name" >> virulence.log

      # Run virulence analysis
      abricate \
        --db $db_to_use \
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
        # Determine risk level
        risk = "Low";
        if ($9 >= 90 && $10 >= 90) risk = "High";
        else if ($9 >= 70 && $10 >= 70) risk = "Medium";

        print sample, $5, $13, $9, $10, risk;
      }' "${output_tsv}.tmp" > "$output_tsv"
      rm "${output_tsv}.tmp"

      # Generate HTML version
      awk -v header="$HTML_HEADER" -v footer="$HTML_FOOTER" '
      BEGIN {
        print header;
      }
      NR==1 { next } # Skip header
      {
        risk_class = "risk-low";
        if ($6 == "High") risk_class = "risk-high";
        else if ($6 == "Medium") risk_class = "risk-medium";

        printf "<tr>";
        printf "<td>%s</td>", $1;  # Sample
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
    done

    # Generate combined report
    echo "Generating combined report..." >> virulence.log

    # Create combined TSV
    echo -e "Sample\tVirulence_Factor\tProduct\t%Coverage\t%Identity\tRisk_Level" > virulence_results/combined_virulence.tsv

    # Process all individual TSV files (skip header line from each)
    for tsv_file in virulence_results/*_virulence.tsv; do
      tail -n +2 "$tsv_file" >> virulence_results/combined_virulence.tsv
    done

    # Generate combined HTML
    {
      echo "$HTML_HEADER"

      # Process combined TSV
      while IFS=$'\t' read -r sample vf product coverage identity risk; do
        # Skip header line
        if [ "$sample" == "Sample" ]; then
          continue
        fi

        # Determine if this is an error row
        if [[ "$vf" == *"ERROR"* ]]; then
          echo "<tr style='background-color: #ffdddd'>"
          echo "<td>$sample</td>"
          echo "<td colspan='5'>Virulence analysis failed</td>"
          echo "</tr>"
        else
          risk_class="risk-low"
          if [ "$risk" == "High" ]; then
            risk_class="risk-high"
          elif [ "$risk" == "Medium" ]; then
            risk_class="risk-medium"
          fi

          echo "<tr>"
          echo "<td>$sample</td>"
          echo "<td class=\"virulence\">$vf</td>"
          echo "<td>$product</td>"
          echo "<td>$coverage</td>"
          echo "<td>$identity</td>"
          echo "<td class=\"$risk_class\">$risk</td>"
          echo "</tr>"
        fi
      done < virulence_results/combined_virulence.tsv

      echo "$HTML_FOOTER"
    } > html_results/combined_virulence.html

    echo "Virulence analysis completed at $(date)" >> virulence.log
    echo "Output files:" >> virulence.log
    ls -lh virulence_results/* html_results/* >> virulence.log 2>&1 || true
  >>>

  runtime {
    docker: "staphb/abricate:latest"
    cpu: cpu
    memory: "4G"
    disks: "local-disk 50 HDD"
    preemptible: 2
    continueOnReturnCode: true
    timeout: "12 hours"
  }

  output {
    Array[File] virulence_reports = glob("virulence_results/*_virulence.tsv")
    Array[File] html_reports = glob("html_results/*_virulence.html")
    File combined_report = "virulence_results/combined_virulence.tsv"
    File combined_html = "html_results/combined_virulence.html"
    File virulence_log = "virulence.log"
  }
}


task BLAST_ANALYSIS {
  input {
    Array[File]+ contig_fastas
    String blast_db
    Boolean use_local_blast = false
    File? local_blast_db
    Int max_target_seqs = 10
    Float evalue = 0.000001
    Int min_contig_length = 200
    Boolean do_blast = true
    Int cpu = 8
    Int memory_gb = 12
    Boolean skip = false
    Boolean debug = false
    Boolean parse_seqids = false
    Int max_retries_per_sample = 1
    Int retry_delay_seconds = 30
  }

  command <<<
    # Always create a status/skip log so the output is never missing
    : > skip_blast.log
    echo "BLAST_ANALYSIS started on $(date)" >> skip_blast.log

    set -euo pipefail

    if [ "~{skip}" == "true" ]; then
      echo "Skipping BLAST analysis as requested" >> skip_blast.log
      mkdir -p blast_results
      : > sample_ids.txt

      for contig_file in ~{sep=' ' contig_fastas}; do
        sample_id=$(basename "$contig_file" | cut -d'.' -f1)
        sample_dir="blast_results/${sample_id}"
        mkdir -p "$sample_dir"

        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle" > "${sample_dir}/blast_results.tsv"
        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle" > "${sample_dir}/top_hits.tsv"
        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle\tspecies" > "${sample_dir}/top_hits_with_species.tsv"

        # Minimal HTML (no big headers) so MERGE_REPORTS won't show duplicate headings
        cat > "${sample_dir}/report.html" <<EOF
<!DOCTYPE html><html><head><meta charset="utf-8">
<style>body{font-family:Arial;margin:0}.note{padding:10px}</style></head>
<body><div class="note">BLAST analysis was skipped for sample ${sample_id}.</div></body></html>
EOF

        echo "$sample_id" >> sample_ids.txt
      done

      echo "Skipping completed at $(date)" >> skip_blast.log
      exit 0
    fi

    export TMPDIR="$(pwd)/blast_tmp"
    mkdir -p "$TMPDIR"
    export PATH="/root/edirect:$PATH"
    DEBUG=~{debug}

    log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2; }
    debug_log() { [ "$DEBUG" = true ] && echo "[DEBUG][$(date '+%H:%M:%S')] $1" >&2; }
    cleanup() { rm -rf "$TMPDIR"; }
    trap cleanup EXIT

    check_ncbi_tools() {
      local tools=("efetch" "esearch" "elink" "xtract")
      local missing=0
      for tool in "${tools[@]}"; do
        command -v "$tool" >/dev/null || { log "WARNING: $tool not found in PATH"; missing=1; }
      done
      [ $missing -eq 0 ]
    }

    if [ "~{use_local_blast}" = "true" ] && [ -n "~{local_blast_db}" ]; then
      log "Using local BLAST database"
      if [ ! -s "~{local_blast_db}" ]; then
        log "Error: Input FASTA for local DB is empty or missing"; exit 1
      fi

      for attempt in {1..3}; do
        makeblastdb -in "~{local_blast_db}" -dbtype nucl -parse_seqids \
                    -title "~{blast_db}" -out "~{blast_db}" \
                    -logfile "makeblastdb.log" 2>&1 && break || {
          log "Attempt $attempt failed to build BLAST DB"
          [ $attempt -lt 3 ] && { sleep 10; rm -f "~{blast_db}".n* "makeblastdb.log"; } || { log "Failed to build local DB"; exit 1; }
        }
      done
      BLAST_DB="~{blast_db}"
    else
      log "Using remote or prebuilt DB: ~{blast_db}"
      BLAST_DB="~{blast_db}"
    fi

    NCBI_TOOLS_AVAILABLE=false
    if [ "~{use_local_blast}" = "true" ] && [ -n "~{local_blast_db}" ]; then
      if check_ncbi_tools; then NCBI_TOOLS_AVAILABLE=true; fi
    fi

    clean_accession() { echo "$1" | sed -E 's/^[^|]*\|([^|]+)\|.*/\1/' | sed 's/\..*//'; }

    get_species() {
      local acc="$1"
      local clean_acc=$(clean_accession "$acc")
      local species="Unknown"
      local header_text=""

      if [ "~{use_local_blast}" = "true" ] && [ -n "~{local_blast_db}" ]; then
        for attempt in {1..3}; do
          header_text=$(blastdbcmd -db "$BLAST_DB" -entry "$clean_acc" -outfmt "%t" 2>&1) && break || sleep 2
        done
        species=$(echo "$header_text" | grep -oP '(\[organism=\K[^]]+|^[^[]+? \K[^[ ]+)' | head -1)
        [ -z "$species" ] && species="Unknown"
        species=$(echo "$species" | sed 's/[]\[,;].*//' | xargs | tr -d '\n')
      fi

      if [ "$species" = "Unknown" ] && [ "$NCBI_TOOLS_AVAILABLE" = true ]; then
        for ((i=1; i<=3; i++)); do
          species=$(efetch -db nuccore -id "$clean_acc" -format docsum 2>/dev/null | \
                   xtract -pattern DocumentSummary -element Organism 2>/dev/null | head -1 | tr -d '\n') && [ -n "$species" ] && break
          sleep 2
        done
      fi
      echo "${species:-Unknown}"
    }

    annotate_blast_results() {
      local input="$1"; local output="$2"
      echo -e "$(head -n1 "$input")\tspecies" > "$output"
      tail -n +2 "$input" | while IFS=$'\t' read -r -a fields; do
        local acc="${fields[1]}"; local species=$(get_species "$acc")
        printf "%s\t%s\n" "$(IFS=$'\t'; echo "${fields[*]}")" "$species" >> "$output"
      done
    }

    generate_html_report() {
      local dir="$1"; local name="$2"
      local report="$dir/report.html"; local input_tsv="$dir/top_hits_with_species.tsv"

      # If no table, emit a tiny message (no headers) so MERGE_REPORTS doesn't duplicate headings
      if [ ! -s "$input_tsv" ]; then
        cat > "$report" <<'EOF'
<!DOCTYPE html><html><head><meta charset="utf-8">
<style>body{font-family:Arial;margin:0}.note{padding:10px}</style></head>
<body><div class="note">No significant BLAST hits.</div></body></html>
EOF
        return 0
      fi

      # Emit only the table; MERGE_REPORTS wraps it with its own section heading
      cat > "$report" <<'EOF'
<!DOCTYPE html><html><head><meta charset="utf-8">
<style>
  body{font-family:Arial;margin:0}
  table{border-collapse:collapse;width:100%}
  th,td{border:1px solid #ddd;padding:8px;text-align:left}
  th{background:#f2f2f2;position:sticky;top:0}
  .species{font-weight:bold}
</style></head><body>
<table><thead><tr>
EOF
      head -n 1 "$input_tsv" | awk -F'\t' '{for(i=1;i<=NF;i++) printf "<th>%s</th>\n",$i}' >> "$report"
      cat >> "$report" <<'EOF'
</tr></thead><tbody>
EOF
      tail -n +2 "$input_tsv" | awk -F'\t' '{
        printf "<tr>";
        for(i=1;i<=NF;i++){
          if(i==NF) printf "<td class=\"species\">%s</td>",$i; else printf "<td>%s</td>",$i
        }
        printf "</tr>\n"
      }' >> "$report"
      echo "</tbody></table></body></html>" >> "$report"
      log "Generated HTML report for $name at $report"
    }

    process_sample() {
      local contig_file="$1"
      local sample_id=$(basename "$contig_file" | cut -d'.' -f1)
      local sample_dir="blast_results/${sample_id}"
      mkdir -p "$sample_dir"
      log "Processing sample: $sample_id"

      local filtered_contig="${sample_dir}/filtered_contigs.fa"
      if [ ! -s "$contig_file" ]; then
        echo ">empty_sequence" > "$filtered_contig"
      else
        awk -v min_len="~{min_contig_length}" 'BEGIN{RS=">";FS="\n"} NR>1{seq=""; for(i=2;i<=NF;i++) seq=seq $i; if(length(seq)>=min_len){print ">"$1; for(i=2;i<=NF;i++) print $i}}' "$contig_file" > "$filtered_contig" || echo ">empty_sequence" > "$filtered_contig"
      fi

      local blast_out="${sample_dir}/blast_results.tsv"
      local blast_log="${sample_dir}/blast.log"
      for attempt in {1..~{max_retries_per_sample}}; do
        blastn -query "$filtered_contig" -db "$BLAST_DB" -task blastn \
               -word_size 28 -reward 1 -penalty -2 -gapopen 2 -gapextend 1 \
               -outfmt "6 std qlen slen stitle" -out "$blast_out" \
               -evalue ~{evalue} -max_target_seqs ~{max_target_seqs} \
               -num_threads ~{cpu} > "$blast_log" 2>&1 && break || {
          log "BLAST attempt $attempt failed for $sample_id"
          [ $attempt -eq ~{max_retries_per_sample} ] && echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle" > "$blast_out" || sleep ~{retry_delay_seconds}
        }
      done

      local top_hits="${sample_dir}/top_hits.tsv"
      { echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle"
        sort -t$'\t' -k12,12g "$blast_out" | head -10; } > "$top_hits" || true

      local annotated_hits="${sample_dir}/top_hits_with_species.tsv"
      annotate_blast_results "$top_hits" "$annotated_hits" || cp "$top_hits" "$annotated_hits"

      generate_html_report "$sample_dir" "$sample_id"
      echo "$sample_id" >> sample_ids.txt
    }

    mkdir -p blast_results
    : > sample_ids.txt
    for contig_file in ~{sep=' ' contig_fastas}; do
      process_sample "$contig_file" || { log "Error on $contig_file; continuing"; continue; }
    done
  >>>

  runtime {
    docker: "gmboowa/blast-analysis:1.9.3"
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk 200 HDD"
    maxRetries: 3
    preemptible: 2
    internet: true
    timeout: "24 hours"
  }

  output {
    Array[File] blast_results = glob("blast_results/*/blast_results.tsv")
    Array[File] blast_top10 = glob("blast_results/*/top_hits.tsv")
    Array[File] blast_annotated = glob("blast_results/*/top_hits_with_species.tsv")
    Array[File] blast_reports = glob("blast_results/*/report.html")
    Array[File] blast_logs = glob("blast_results/*/blast.log")
    Array[File] system_logs = glob("blast_results/*/makeblastdb.log")
    Array[String] sample_ids = read_lines("sample_ids.txt")
    File skip_log = "skip_blast.log"
  }
}

task TREE_VISUALIZATION {
  input {
    File input_tree
    Int width = 1200
    Int height = 1000
    String image_format = "png"
    Int font_size = 10
    Int title_font_size = 12
    String layout = "circular"  # Options: circular, rectangular, or fan
    Float branch_thickness = 1.5
    Boolean show_branch_lengths = false
    Boolean show_scale = true
    String color_scheme = "standard"  # Options: standard, gradient, or categorical
    String? tree_title
  }

  command <<<
    set -euo pipefail

    # Create output directory
    mkdir -p output

    # Get tree type from filename
    tree_basename=$(basename "~{input_tree}" .nwk)
    tree_type=""
    color_scheme_val="~{color_scheme}"

    if [[ "$tree_basename" == *"accessory"* ]]; then
      tree_type="accessory"
      # If no explicit color scheme was passed, default to categorical for accessory
      if [[ -z "${color_scheme_val}" || "${color_scheme_val}" == "standard" ]]; then
        color_scheme_val="categorical"
      fi
    elif [[ "$tree_basename" == *"core"* ]]; then
      tree_type="core"
      # If no explicit color scheme was passed, default to gradient for core
      if [[ -z "${color_scheme_val}" || "${color_scheme_val}" == "standard" ]]; then
        color_scheme_val="gradient"
      fi
    fi

    # Title: prefer user-provided tree_title, otherwise build a sensible default in bash
    if [[ -n "~{tree_title}" ]]; then
      TITLE="~{tree_title}"
    else
      if [[ "$tree_type" == "accessory" ]]; then
        TITLE="Accessory Genes Phylogenetic Tree"
      else
        TITLE="Core Genes Phylogenetic Tree"
      fi
    fi

    export TITLE
    export color_scheme_val

    python3 <<EOF
  import os
  from ete3 import Tree, TreeStyle, NodeStyle, TextFace
  import colorsys
  import random

  # Inputs from bash/WDL
  input_tree = "~{input_tree}"
  width = ~{width}
  height = ~{height}
  show_branch_lengths = ~{show_branch_lengths}
  show_scale = ~{show_scale}
  layout = "~{layout}".lower()
  branch_thickness = ~{branch_thickness}
  title_text = os.environ.get("TITLE", "Phylogenetic Tree")
  color_scheme = os.environ.get("color_scheme_val", "standard").lower()

  # Load the tree
  tree = Tree(input_tree)

  # TreeStyle
  ts = TreeStyle()
  ts.show_leaf_name = True
  ts.show_branch_length = show_branch_lengths
  ts.show_scale = show_scale
  ts.branch_vertical_margin = 10
  ts.scale = ~{width / 6}  # use / not // in WDL
  ts.title.add_face(TextFace(title_text, fsize=~{title_font_size}, fgcolor="black"), column=0)

  # Layout
  if layout == "circular":
      ts.mode = "c"
      ts.arc_start = -180
      ts.arc_span = 240
      ts.rotation = 90
  elif layout == "fan":
      ts.mode = "c"
      ts.arc_start = 0
      ts.arc_span = 360
  else:
      ts.mode = "r"
      ts.root_opening_factor = 1.0

  # Color schemes
  if color_scheme == "gradient":
      # Gradient color from blue to red based on branch length
      max_dist = max((node.dist for node in tree.traverse()), default=1.0)
      for node in tree.traverse():
          if node.dist > 0:
              ratio = min(node.dist / max_dist, 1.0)
              r, g, b = colorsys.hsv_to_rgb(0.66 * (1 - ratio), 0.7, 0.8)
              node.img_style["hz_line_color"] = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
  elif color_scheme == "categorical":
      random.seed(42)
      for node in tree.traverse():
          if not node.is_leaf():
              hue = random.uniform(0, 0.8)
              r, g, b = colorsys.hsv_to_rgb(hue, 0.8, 0.8)
              node.img_style["hz_line_color"] = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"

  # Common node style
  for n in tree.traverse():
      n.img_style["size"] = 0
      n.img_style["hz_line_width"] = branch_thickness
      n.img_style["vt_line_width"] = branch_thickness
      if n.is_leaf():
          n.img_style["fgcolor"] = "black"
          n.img_style["shape"] = "square"
          n.img_style["size"] = 5

  # Leaf font size
  leaf_style = NodeStyle()
  leaf_style["size"] = ~{font_size}
  for leaf in tree.iter_leaves():
      leaf.set_style(leaf_style)

  # Render
  output_file = f"output/{os.path.basename(input_tree)[:-4]}.~{image_format}"
  tree.render(output_file, w=width, h=height, units="px", tree_style=ts, dpi=300)

  # Log
  with open(f"output/{os.path.basename(input_tree)[:-4]}.log", "w") as logfile:
      logfile.write("Tree visualization parameters:\\n")
      logfile.write(f"Input tree: {input_tree}\\n")
      logfile.write(f"Dimensions: {width}x{height}px\\n")
      logfile.write(f"Layout: {layout}\\n")
      logfile.write(f"Color scheme: {color_scheme}\\n")
      logfile.write(f"Branch thickness: {branch_thickness}\\n")
      logfile.write(f"Output image: {output_file}\\n")
  EOF
  >>>


  runtime {
    docker: "gmboowa/ete3-render:1.14"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
  }

  output {
    File final_image = "output/~{basename(input_tree, '.nwk')}.~{image_format}"
    File render_log = "output/~{basename(input_tree, '.nwk')}.log"
  }
}
task MERGE_REPORTS {
  input {
    # Section HTMLs (optional)
    File? trimming_report_html
    File? fastqc_summary_html
    File? assembly_stats_html
    File? annotation_summary_html
    File? pangenome_report_html
    File? gene_heatmap_png
    File? mlst_combined_html
    File? variant_summary_html
    File? amr_combined_html
    File? mge_combined_html
    File? virulence_combined_html

    # Direct-file inputs
    File? pangenome_report
    File? gene_heatmap

    # Collections
    Array[File] blast_reports = []
    Array[File] tree_images = []

    # Meta
    String workflow_name = "rMAP Analysis Pipeline"
    String version = "2.0"
    String footer_sentence = ""
    Boolean skip = false
  }

  command <<<
    #!/usr/bin/env bash
    set -euo pipefail

    # Always ensure this file exists so Cromwell can map it
    : > skip_report.log

    # Always start clean
    rm -rf final_report default_report.tgz final_report.tgz || true
    mkdir -p final_report/assets/images

    RUN_DATE="$(date)"
    RUN_HEADER="This run was generated using rMAP version ~{version} performed on ${RUN_DATE}… Thank you for using this pipeline!!"

    if [ "~{skip}" == "true" ]; then
      echo "Skipping report generation as requested" > skip_report.log

      cat > final_report/index.html <<'EOF'
<!DOCTYPE html>
<html lang="en"><head>
  <meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Report Generation Skipped</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; }
    .skip-notice { background:#fff3cd;border:1px solid #ffeeba;padding:20px;border-radius:5px;margin:20px 0;text-align:center; }
  </style>
</head>
<body>
  <h1>Report Generation Skipped</h1>
  <div class="skip-notice">
    <p>This report was not generated because the workflow was configured to skip this step.</p>
    <p>Generated: DATE_PLACEHOLDER</p>
  </div>
</body>
</html>
EOF
      # Stamp date
      sed -i "s/DATE_PLACEHOLDER/${RUN_DATE}/" final_report/index.html

      tar -czf final_report.tgz final_report
      exit 0
    fi

    # Prepare a tiny default archive in case something fails later
    echo "<html><body><h1>Report Generation Failed</h1></body></html>" > final_report/index.html
    tar -czf default_report.tgz final_report/index.html

    # Helper to copy if present
    cp_if_exists() {
      local src="$1"
      local dst_dir="$2"
      if [ -n "${src}" ] && [ -f "${src}" ]; then
        cp "${src}" "${dst_dir}/"
      fi
    }

    # Heatmap image (two possible inputs)
    GH_BASENAME=""
    if [ -n "~{gene_heatmap}" ] && [ -f "~{gene_heatmap}" ]; then
      cp "~{gene_heatmap}" final_report/assets/images/
      GH_BASENAME="$(basename "~{gene_heatmap}")"
    elif [ -n "~{gene_heatmap_png}" ] && [ -f "~{gene_heatmap_png}" ]; then
      cp "~{gene_heatmap_png}" final_report/assets/images/
      GH_BASENAME="$(basename "~{gene_heatmap_png}")"
    fi

    # Copy tree images (if any)
    for img in ~{sep=' ' tree_images}; do
      if [ -f "$img" ]; then
        cp "$img" final_report/assets/images/
      fi
    done

    # HTML embedding function (safe: checks file existence)
    embed_html() {
      local label="$1"
      local file="$2"
      local out="$3"

      {
        echo "          <div class=\"report-card\">"
        echo "            <h3>${label}</h3>"
        echo "            <div class=\"report-content\">"
        if [ -n "${file}" ] && [ -f "${file}" ]; then
          if grep -q "<table" "$file" 2>/dev/null; then
            sed -n '/<table[ >]/,/<\/table>/p' "$file"
          else
            cat "$file"
          fi
        else
          echo "<p><em>No ${label} available.</em></p>"
        fi
        echo "            </div>"
        echo "          </div>"
      } >> "$out"
    }

    # --------- Build the main HTML ----------
    cat > final_report/index.html <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>~{workflow_name} - Comprehensive Report</title>
  <style>
    :root { --primary:#3498db; --secondary:#2c3e50; --light:#f8f9fa; --text:#333; --accent:#fff3cd; }
    body { font-family:Segoe UI,Tahoma,Arial,sans-serif; margin:0; background:#f5f5f5; color:var(--text) }
    .container { max-width:1200px; margin:0 auto; padding:20px }
    .header { background:var(--secondary); color:#fff; padding:24px; border-radius:6px; margin-bottom:24px }
    h1 { margin:0; font-size:2rem }
    .meta { opacity:.9 }
    .run-banner { margin-top:10px; background:var(--accent); color:#333; padding:10px 12px; border-radius:6px; font-weight:600; border:1px solid #e6cf8b }
    .layout { display:flex; gap:24px }
    .nav { position:sticky; top:16px; width:260px; background:var(--light); padding:16px; border-radius:6px }
    .nav a { display:block; padding:8px 10px; border-radius:4px; text-decoration:none; color:var(--secondary) }
    .nav a:hover { background:var(--primary); color:#fff }
    .content { flex:1; background:#fff; padding:20px; border-radius:6px }
    .section { margin-bottom:36px; border-bottom:1px solid #eee; padding-bottom:20px }
    .section:last-child { border:none }
    .report-card { background:var(--light); padding:16px; border-radius:6px; margin:12px 0 }
    .image-center { text-align:center; margin:12px 0 }
    img { max-width:100%; height:auto; border:1px solid #ddd; border-radius:4px }
    .tabs { display:flex; gap:8px; margin:12px 0 }
    .tab-btn { border:0; padding:8px 12px; border-radius:4px; background:#eee; cursor:pointer }
    .tab-btn.active { background:var(--primary); color:#fff }
    .tab { display:none }
    .tab.active { display:block }
  </style>
</head>
<body>
<div class="container">
  <div class="header">
    <h1>~{workflow_name}</h1>
    <div class="meta">Version ~{version} • Generated <span id="gen-date"></span></div>
    <div class="run-banner">${RUN_HEADER}</div>
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
          <p>This report aggregates results from rMAP modules (QC, assembly/annotation, pangenome/phylogeny, AMR/MGE/virulence, and BLAST).</p>
        </div>
      </div>
EOF

    # QC
    cat >> final_report/index.html <<'EOF'
      <div id="qc" class="section">
        <h2>Quality Control</h2>
EOF
    embed_html "Read trimming" "~{trimming_report_html}" "final_report/index.html"
    embed_html "FastQC / MultiQC summary" "~{fastqc_summary_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # Assembly
    cat >> final_report/index.html <<'EOF'
      <div id="assembly" class="section">
        <h2>Assembly</h2>
EOF
    embed_html "Assembly statistics" "~{assembly_stats_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # Annotation
    cat >> final_report/index.html <<'EOF'
      <div id="annotation" class="section">
        <h2>Annotation</h2>
EOF
    embed_html "Annotation summary" "~{annotation_summary_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # Pangenome
    cat >> final_report/index.html <<'EOF'
      <div id="pangenome" class="section">
        <h2>Pangenome</h2>
EOF
    if [ -n "~{pangenome_report}" ] && [ -f "~{pangenome_report}" ]; then
      embed_html "Pangenome report" "~{pangenome_report}" "final_report/index.html"
    else
      embed_html "Pangenome report" "~{pangenome_report_html}" "final_report/index.html"
    fi

    if [ -n "${GH_BASENAME:-}" ]; then
      cat >> final_report/index.html <<EOF
        <div class="report-card">
          <h3>Gene Presence/Absence Heatmap</h3>
          <div class="image-center">
            <img src="assets/images/${GH_BASENAME}" alt="Gene presence/absence heatmap">
          </div>
        </div>
EOF
    fi
    echo "      </div>" >> final_report/index.html

    # MLST
    cat >> final_report/index.html <<'EOF'
      <div id="mlst" class="section">
        <h2>MLST</h2>
EOF
    embed_html "Combined MLST" "~{mlst_combined_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # Variants
    cat >> final_report/index.html <<'EOF'
      <div id="variants" class="section">
        <h2>Variants</h2>
EOF
    embed_html "Variant summary" "~{variant_summary_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # AMR
    cat >> final_report/index.html <<'EOF'
      <div id="amr" class="section">
        <h2>Antimicrobial Resistance</h2>
EOF
    embed_html "Combined AMR" "~{amr_combined_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # MGE
    cat >> final_report/index.html <<'EOF'
      <div id="mge" class="section">
        <h2>Mobile Genetic Elements</h2>
EOF
    embed_html "Combined MGE" "~{mge_combined_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # Virulence
    cat >> final_report/index.html <<'EOF'
      <div id="virulence" class="section">
        <h2>Virulence</h2>
EOF
    embed_html "Combined virulence" "~{virulence_combined_html}" "final_report/index.html"
    echo "      </div>" >> final_report/index.html

    # BLAST tabs (only for existing files)
    cat >> final_report/index.html <<'EOF'
      <div id="blast" class="section">
        <h2>BLAST</h2>
        <div class="tabs">
EOF
    i=0
    for br in ~{sep=' ' blast_reports}; do
      if [ -f "$br" ]; then
        sname=$(basename "$(dirname "$br")")
        cls=""
        if [ $i -eq 0 ]; then cls="active"; fi
        echo "<button class=\"tab-btn ${cls}\" onclick=\"openTab(event,'tab${i}')\">${sname}</button>" >> final_report/index.html
        i=$((i+1))
      fi
    done
    echo "        </div>" >> final_report/index.html

    i=0
    for br in ~{sep=' ' blast_reports}; do
      if [ -f "$br" ]; then
        cls=""
        if [ $i -eq 0 ]; then cls="active"; fi
        echo "<div id=\"tab${i}\" class=\"tab ${cls}\">" >> final_report/index.html
        if grep -q "<table" "$br" 2>/dev/null; then
          sed -n '/<table[ >]/,/<\/table>/p' "$br" >> final_report/index.html
        else
          echo "<p><em>BLAST report available: $(basename "$br")</em></p>" >> final_report/index.html
        fi
        echo "</div>" >> final_report/index.html
        i=$((i+1))
      fi
    done
    echo "      </div>" >> final_report/index.html

    # Phylogeny
    cat >> final_report/index.html <<'EOF'
      <div id="phylogeny" class="section">
        <h2>Phylogeny</h2>
EOF
    for img in ~{sep=' ' tree_images}; do
      if [ -f "$img" ]; then
        b=$(basename "$img")
        cat >> final_report/index.html <<EOF
        <div class="report-card">
          <h3>$b</h3>
          <div class="image-center"><img src="assets/images/$b" alt="Phylogenetic tree"></div>
        </div>
EOF
      fi
    done
    echo "      </div>" >> final_report/index.html

    # Footer + closing (footer text stamped below)
    cat >> final_report/index.html <<'EOF'
    </div>
  </div>
  <div class="footer-banner" style="text-align:center;color:#333;margin:24px 0;padding:10px 12px;background:#fff3cd;border:1px solid #e6cf8b;border-radius:6px;font-weight:600">
    FOOTER_SENTENCE_PLACEHOLDER
  </div>
</div>
<script>
  document.getElementById('gen-date').textContent = new Date().toLocaleString();
  function openTab(evt, id){
    document.querySelectorAll('.tab').forEach(t=>t.classList.remove('active'));
    document.querySelectorAll('.tab-btn').forEach(b=>b.classList.remove('active'));
    document.getElementById(id).classList.add('active');
    evt.currentTarget.classList.add('active');
  }
</script>
</body>
</html>
EOF

    # Build footer sentence with dynamic date/version, honoring user text if provided.
    FS="~{footer_sentence}"
    if [ -z "${FS// }" ]; then
      FS="This run was generated using rMAP version ~{version} performed on ${RUN_DATE}… Thank you for using this pipeline!!"
    else
      FS="${FS//{RUN_DATE}/${RUN_DATE}}"
      FS="${FS//{VERSION}/~{version}}"
      FS="${FS//2021-06-14/${RUN_DATE}}"
    fi
    # Stamp footer into HTML
    sed -i "s|FOOTER_SENTENCE_PLACEHOLDER|${FS}|" final_report/index.html

    # Create archive; if tar fails, fall back to default
    tar -czf final_report.tgz final_report || cp default_report.tgz final_report.tgz
  >>>

  runtime {
    docker: "ubuntu:20.04"
    cpu: 2
    memory: "2 GB"
    disks: "local-disk 50 HDD"
    continueOnReturnCode: true
    maxRetries: 1
  }

  output {
    File final_report_html = "final_report/index.html"
    File final_report_tgz = "final_report.tgz"
    File skip_log = "skip_report.log"
  }
}
