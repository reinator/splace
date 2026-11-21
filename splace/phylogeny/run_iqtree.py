import logging
import os
import re
import shutil
import subprocess
from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def check_iqtree():
    """Return the IQ-TREE executable on PATH: iqtree3 (IQ-TREE 3, installed by
    environment.yml), then iqtree2 or iqtree; None when there is none."""
    for name in ("iqtree3", "iqtree2", "iqtree"):
        if shutil.which(name):
            return name
    return None

def _format_taxon_label(raw_name):
    """Parse 'Genus_species_AccessionID_UID' into a display label with italic species.

    Uses matplotlib mathtext for italic rendering. Spaces must be escaped
    as '\\ ' inside math mode, otherwise they are ignored.

    Expected format: Genus_species_NC_123456_00001 or Genus_species_SomeID_00001
    The last 5-digit segment is always the unique record ID.
    """
    # Extract trailing 5-digit UID
    uid_match = re.match(r"^(.+)_(\d{5})$", raw_name)
    if uid_match:
        base = uid_match.group(1)
        uid = uid_match.group(2)
    else:
        base = raw_name
        uid = ""

    # Try to extract accession from base: Genus_species_NC_123456 or Genus_species_MN123456.1
    acc_match = re.match(r"^(.+?)_(([A-Z]{1,2}_?\d{5,}(?:\.\d+)?))$", base)
    if acc_match:
        species_part = acc_match.group(1).replace("_", r"\ ")
        accession = acc_match.group(2)
        label = f"$\\it{{{species_part}}}$ ({accession})"
    else:
        species_part = base.replace("_", r"\ ")
        label = f"$\\it{{{species_part}}}$"

    if uid:
        label += f" [{uid}]"
    return label

def _build_presence_report(all_taxa, gene_data, output_dir):
    """Build a gene presence/absence text report and heatmap."""
    sorted_taxa = sorted(all_taxa)
    gene_names = [g[0] for g in gene_data]
    num_genes = len(gene_names)
    num_taxa = len(sorted_taxa)

    # Build presence matrix
    matrix = []
    for taxon in sorted_taxa:
        row = [1 if taxon in records else 0 for _, _, records in gene_data]
        matrix.append(row)

    # --- Text report (console + file) ---
    taxon_col_w = max(len(t) for t in sorted_taxa + ["Taxon", "Present"]) + 2
    gene_col_w = {g: max(len(g), 3) + 2 for g in gene_names}
    total_col_w = max(len(f"0/{num_genes}"), len("Total")) + 2

    header = f"{'Taxon':<{taxon_col_w}}"
    for g in gene_names:
        header += f" | {g:^{gene_col_w[g]}}"
    header += f" | {'Total':^{total_col_w}}"
    sep = "-" * len(header)

    lines = ["", "Gene Presence/Absence Report", "=" * len(header), header, sep]

    gene_totals = [0] * num_genes
    for idx, taxon in enumerate(sorted_taxa):
        row_str = f"{taxon:<{taxon_col_w}}"
        count = 0
        for i, g in enumerate(gene_names):
            present = matrix[idx][i]
            symbol = "+" if present else "-"
            row_str += f" | {symbol:^{gene_col_w[g]}}"
            if present:
                count += 1
                gene_totals[i] += 1
        row_str += f" | {count}/{num_genes}".rjust(total_col_w + 3)
        lines.append(row_str)

    lines.append(sep)
    summary = f"{'Present':<{taxon_col_w}}"
    for i, g in enumerate(gene_names):
        val = f"{gene_totals[i]}/{num_taxa}"
        summary += f" | {val:^{gene_col_w[g]}}"
    summary += f" |"
    lines.extend([summary, ""])

    table_text = "\n".join(lines)
    logging.info(table_text)

    report_path = os.path.join(output_dir, "gene_presence_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(table_text + "\n")
    logging.info(f"Gene presence report saved to: {report_path}")

    # --- Heatmap ---
    df = pd.DataFrame(matrix, index=sorted_taxa, columns=gene_names)

    # Format taxon labels with italic species names (UID already embedded in taxon name)
    formatted_labels = [_format_taxon_label(t) for t in sorted_taxa]

    # Dynamic figure sizing
    cell_w, cell_h = 0.5, 0.35
    fig_w = max(8, num_genes * cell_w + 4)
    fig_h = max(4, num_taxa * cell_h + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = sns.color_palette(["#f26c41", "#4676b4"])
    sns.heatmap(
        df,
        ax=ax,
        cmap=cmap,
        cbar=False,
        linewidths=0.5,
        linecolor="white",
        square=True,
        xticklabels=True,
        yticklabels=formatted_labels,
    )

    ax.set_xlabel("Genes", fontsize=12, fontweight="bold")
    ax.set_ylabel("")
    ax.set_title("Gene Presence / Absence", fontsize=14, fontweight="bold", pad=15)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    plt.xticks(rotation=45, ha="left", fontsize=9)
    plt.yticks(fontsize=9)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4676b4", edgecolor="white", label="Present"),
        Patch(facecolor="#f26c41", edgecolor="white", label="Absent"),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
        fontsize=10,
    )

    plt.tight_layout()
    heatmap_path = os.path.join(output_dir, "gene_presence_heatmap.png")
    fig.savefig(heatmap_path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Gene presence heatmap saved to: {heatmap_path}")

def prepare_supermatrix(trimmed_files, output_dir, allow_missing=False, report_dir=None):
    """
    Concatenate trimmed alignments into a supermatrix (NEXUS format) with partition block.

    Args:
        trimmed_files (list): List of paths to trimmed FASTA files.
        output_dir (str): Directory to save the supermatrix.
        allow_missing (bool): If True, fill missing genes with '?' characters.
            If False, remove genes that are not present in all taxa.

    Returns:
        str: Path to the generated supermatrix file.
    """
    logging.info("Preparing supermatrix...")

    # Sort files to ensure deterministic order
    trimmed_files.sort()

    all_taxa = set()
    gene_data = [] # List of tuples: (gene_name, length, sequence_dict)

    # 1. Read all files
    for fasta_path in trimmed_files:
        gene_name = os.path.basename(fasta_path).replace("_aligned_trimmed.fasta", "").replace("_trimmed.fasta", "").replace(".fasta", "")

        records = {}
        length = 0
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                if length == 0:
                    length = len(record.seq)
                elif len(record.seq) != length:
                    logging.warning(f"Sequence length mismatch in {gene_name}. Expected {length}, got {len(record.seq)} for {record.description}")

                records[record.description] = str(record.seq)
                all_taxa.add(record.description)

            if length > 0:
                gene_data.append((gene_name, length, records))
            else:
                logging.warning(f"Gene {gene_name} is empty or invalid. Skipping.")

        except Exception as e:
            logging.error(f"Error reading {fasta_path}: {e}")

    if not gene_data:
        logging.error("No valid data found to create supermatrix.")
        return None

    # Report gene presence/absence matrix
    _build_presence_report(all_taxa, gene_data, report_dir or output_dir)

    # 2. Filter genes based on allow_missing
    if not allow_missing:
        complete_genes = []
        for gene_name, length, records in gene_data:
            missing_taxa = all_taxa - set(records.keys())
            if missing_taxa:
                logging.info(f"Removing gene '{gene_name}': missing from {len(missing_taxa)} taxa ({', '.join(sorted(missing_taxa)[:3])}{'...' if len(missing_taxa) > 3 else ''})")
            else:
                complete_genes.append((gene_name, length, records))

        if not complete_genes:
            logging.error("No genes are present in all taxa. Use --allow-missing to include partial data.")
            return None

        removed_count = len(gene_data) - len(complete_genes)
        if removed_count > 0:
            logging.info(f"Kept {len(complete_genes)}/{len(gene_data)} genes (removed {removed_count} incomplete genes)")
        gene_data = complete_genes

    # 3. Build Supermatrix
    output_nexus = os.path.join(output_dir, "supermatrix.nex")
    total_chars = sum(g[1] for g in gene_data)
    num_taxa = len(all_taxa)

    # When genes were filtered, some taxa may no longer have any data
    # Recompute actual taxa from remaining genes
    if not allow_missing:
        actual_taxa = set()
        for _, _, records in gene_data:
            actual_taxa.update(records.keys())
        all_taxa = actual_taxa
        num_taxa = len(all_taxa)

    try:
        with open(output_nexus, "w", encoding="utf-8") as out:
            out.write("#NEXUS\n")
            out.write("BEGIN DATA;\n")
            out.write(f"DIMENSIONS NTAX={num_taxa} NCHAR={total_chars};\n")
            out.write("FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
            out.write("MATRIX\n")

            sorted_taxa = sorted(list(all_taxa))
            max_name_len = max(len(t) for t in sorted_taxa) + 2

            missing_char = "?" if allow_missing else "-"

            for taxon in sorted_taxa:
                safe_taxon_name = f"'{taxon}'"
                out.write(f"{safe_taxon_name.ljust(max_name_len)} ")

                for _, length, records in gene_data:
                    seq = records.get(taxon, missing_char * length)
                    out.write(seq)
                out.write("\n")

            out.write(";\nEND;\n\n")

            # 4. Write Partition Block (Sets)
            out.write("BEGIN SETS;\n")
            current_pos = 1
            for gene_name, length, _ in gene_data:
                end_pos = current_pos + length - 1
                out.write(f"    CHARSET {gene_name} = {current_pos}-{end_pos};\n")
                current_pos = end_pos + 1
            out.write("END;\n")

        logging.info(f"Supermatrix created at {output_nexus} (Taxa: {num_taxa}, Sites: {total_chars}, Genes: {len(gene_data)}, Missing data: {'allowed' if allow_missing else 'not allowed'})")
        return output_nexus

    except Exception as e:
        logging.error(f"Failed to create supermatrix: {e}")
        return None

def run_iqtree_analysis(supermatrix_path, output_dir, threads=1, bootstrap=1000, model="MFP", extra_args=""):
    """
    Run IQ-TREE analysis on the supermatrix.

    Args:
        supermatrix_path (str): Path to the supermatrix file.
        output_dir (str): Directory for output files.
        threads (int): Number of threads.
        bootstrap (int): Number of ultrafast bootstrap replicates.
        model (str): Substitution model (e.g., "MFP" for ModelFinder).
        extra_args (str): Additional IQ-TREE arguments.

    Returns:
        str: Path to the main tree file (.treefile) or None.
    """
    iqtree_cmd = check_iqtree()
    if not iqtree_cmd:
        logging.error("IQ-TREE not found in PATH. Please install IQ-TREE 3 (conda install -c bioconda iqtree).")
        return None # Critical failure for this step

    prefix = os.path.join(output_dir, "splace_tree")
    
    # Command: iqtree3 -s file.nex -T 8 -B 1000 -m MFP --prefix ... -redo
    # -s: Input (Nexus file determines alignment + partitions)
    # -T: number of threads (IQ-TREE 2 and 3; -nt is the old name)
    # -B: Ultrafast Bootstrap (1000)
    # -m: ModelFinder (MFP)

    cmd = [
        iqtree_cmd,
        "-s", supermatrix_path,
        "-T", str(threads),
        "-B", str(bootstrap),
        "-m", model,
        "--prefix", prefix,
        "-redo" # Overwrite existing
    ]
    if extra_args and extra_args.strip():
        cmd.extend(extra_args.strip().split())
    
    logging.info(f"Starting IQ-TREE analysis with {threads} threads...")
    logging.info(f"Command: {' '.join(cmd)}")
    
    try:
        # IQ-TREE writes to stdout/stderr. Capturing and logging might be noisy, 
        # but useful.
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        if process.returncode == 0:
            logging.info("IQ-TREE analysis completed successfully.")
            tree_file = f"{prefix}.treefile"
            if os.path.exists(tree_file):
                return tree_file
            else:
                logging.warning(f"IQ-TREE finished 0 but {tree_file} not found.")
                return None
        else:
            logging.error(f"IQ-TREE failed with return code {process.returncode}")
            # Log the tail of stderr/stdout
            logging.error(f"IQ-TREE stderr: {process.stderr[-1000:]}")
            return None
            
    except Exception as e:
        logging.error(f"Error executing IQ-TREE: {e}")
        return None

def run_phylogeny_pipeline(trimmed_files, output_dir, threads=1, allow_missing=False, bootstrap=1000, model="MFP", extra_args="", report_dir=None):
    """
    Orchestrate the phylogeny step: Supermatrix -> IQ-TREE
    """
    if not trimmed_files:
        logging.warning("No input files for phylogeny.")
        return None

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True, mode=0o755)

    # Step 1: Supermatrix
    supermatrix = prepare_supermatrix(trimmed_files, output_dir, allow_missing=allow_missing, report_dir=report_dir)
    if not supermatrix:
        return None

    # Step 2: IQ-TREE
    tree_file = run_iqtree_analysis(supermatrix, output_dir, threads, bootstrap=bootstrap, model=model, extra_args=extra_args)

    return tree_file
