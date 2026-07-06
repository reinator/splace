#!/usr/bin/env python3

import asyncio
import argparse
import logging
import re
import shutil
import subprocess
import os
import sys
import yaml

from splace.read_genbank.genbank_handler import FASTA_HEADER_FIELDS

__authors__     = "Renato Oliveira and Luan Rabelo"
__license__     = "GPL-3.0"
__version__     = "4.0.0"
__maintainers__ = "Renato Oliveira and Luan Rabelo"
__email__       = "luan.rabelo@pq.itv.org"
__date__        = "2025/11/15"
__github__      = "itvgenomics/splace"
__status__      = "Stable"
__tool__        = "SPLACE"


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y/%m/%d - %H:%M:%S'
    )


DEFAULT_FASTA_HEADER_CONFIG = {
    "template": "{accession}_{family}_{genus}_{species}",
    "missing_value": "?",
}


def ensure_default_fasta_header_config(config_path):
    if os.path.exists(config_path):
        return

    config_dir = os.path.dirname(config_path)
    
    if config_dir:
        os.makedirs(config_dir, exist_ok=True)

    with open(config_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(DEFAULT_FASTA_HEADER_CONFIG, handle, sort_keys=False)

    logging.info(f"Created default FASTA header configuration at: {config_path}")


def load_fasta_header_config(config_path):
    ensure_default_fasta_header_config(config_path)

    with open(config_path, "r", encoding="utf-8") as handle:
        raw_config = yaml.safe_load(handle) or {}

    if isinstance(raw_config, str):
        config = dict(DEFAULT_FASTA_HEADER_CONFIG)
        config["template"] = raw_config
    elif isinstance(raw_config, dict):
        config = dict(DEFAULT_FASTA_HEADER_CONFIG)
        config.update({key: value for key, value in raw_config.items() if value is not None})
    else:
        raise ValueError("FASTA header configuration must be a YAML string or mapping.")

    template = str(config.get("template", "")).strip()
    missing_value = str(config.get("missing_value", "?") or "?")
    placeholders = set(re.findall(r"\{([a-zA-Z_][a-zA-Z0-9_]*)\}", template))

    if not template:
        raise ValueError("FASTA header template cannot be empty.")
    if not placeholders:
        raise ValueError("FASTA header template must contain at least one placeholder like {accession}.")

    invalid_placeholders = sorted(placeholders - FASTA_HEADER_FIELDS)
    if invalid_placeholders:
        valid_fields = ", ".join(sorted(FASTA_HEADER_FIELDS))
        invalid_fields = ", ".join(invalid_placeholders)
        raise ValueError(
            f"Unsupported FASTA header placeholders: {invalid_fields}. Available fields: {valid_fields}."
        )

    return template, missing_value


def collect_input_files(input_directories):
    genbank_files = set()
    fasta_files = set()

    for input_dir in input_directories:
        if not input_dir:
            continue
        if not os.path.exists(input_dir):
            raise FileNotFoundError(f"Input directory '{input_dir}' not found.")

        for root, _, files in os.walk(input_dir):
            for file_name in files:
                file_path = os.path.join(root, file_name)
                lowered = file_name.lower()
                if lowered.endswith((".gb", ".gbk", ".genbank")):
                    genbank_files.add(file_path)
                elif lowered.endswith((".fasta", ".fa", ".fas")):
                    fasta_files.add(file_path)

    return sorted(genbank_files), sorted(fasta_files)

parser = argparse.ArgumentParser(
    prog=__tool__,
    usage='%(prog)s [options]',
    #description=pipenote_ascii_art,
    epilog=f"In Development...",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=True,
    )

input_group = parser.add_argument_group(
    title="Input Files",
    description="Specify the directory that contains the sequence files to be processed."
)
input_group.add_argument(
    "-i", "--input_dir",
    required=False,
    type=str,
    help="Directory with FASTA (.fasta, .fa, .fas) or GenBank (.gb, .gbk) files. Optional when --ncbi-search-term is used.",
    default=None
)

output_group = parser.add_argument_group(
    title="Output Directory",
    description="Specify the directory where output files will be saved."
)
output_group.add_argument(
    "-o", "--output_dir",
    required=False,
    type=str,
    help="Directory to save output files. Default: 'splace_output/'.",
    default="splace_output/"
)

threads = parser.add_argument_group(
    title="Threads",
    description="Number of threads to use."
    )
threads.add_argument(
    "-t", "--threads",
    default=8,
    type=int,
    help="Number of threads to use. Default is 8."
    )

genbank_group = parser.add_argument_group(
    title="GenBank Options",
    description="Options for processing GenBank files."
)
genbank_group.add_argument(
    "--gb-type",
    choices=["mt", "cp"],
    help="Extract mitochondrial or chloroplast sequences from GenBank files.",
    default="mt"
)
genbank_group.add_argument(
    "--gbif",
    action="store_true",
    help="Query GBIF taxonomy for valid binomial species names found in GenBank records and store the retrieved ranks in the run metadata."
)
genbank_group.add_argument(
    "--fasta-header-config",
    type=str,
    default="fasta_header.yaml",
    help="Path to the YAML file that defines the FASTA header template. Default: 'fasta_header.yaml'."
)

remote_group = parser.add_argument_group(
    title="Remote Retrieval",
    description="Download GenBank genomes from NCBI using a taxonomic name/rank plus genome scope flags."
)
remote_group.add_argument(
    "--ncbi-search-term",
    type=str,
    default=None,
    help="Taxonomic name or rank used to build the NCBI genome query automatically (for example: 'Bufonidae' or 'Coffea')."
)
remote_group.add_argument(
    "--apis-env",
    type=str,
    default="apis.env",
    help="Path to an optional apis.env file containing NCBI_API_KEY and NCBI_EMAIL. Default: 'apis.env'."
)
remote_group.add_argument(
    "--ncbi-download-dir",
    type=str,
    default=None,
    help="Directory where GenBank files downloaded from NCBI will be saved."
)
remote_group.add_argument(
    "--ncbi-complete",
    action="store_true",
    help="Include complete genomes in the automatically generated NCBI query."
)
remote_group.add_argument(
    "--ncbi-partial",
    action="store_true",
    help="Include partial or incomplete genomes in the automatically generated NCBI query."
)
remote_group.add_argument(
    "--ncbi-refseq-only",
    action="store_true",
    help="Restrict the automatically generated NCBI query to RefSeq genomes only."
)
remote_group.add_argument(
    "--download-only",
    action="store_true",
    help="Run only the NCBI search/download stage and stop before marker extraction or downstream analyses."
)

gene_filter_group = parser.add_argument_group(
    title="Gene Filtering",
    description="Filter sequences either by specific gene names or by GenBank feature types."
)
gene_filter_options = gene_filter_group.add_mutually_exclusive_group()
gene_filter_options.add_argument(
    "--genes",
    type=str,
    default=None,
    help="Comma-separated list of gene names (e.g., '12S,16S,COI') or path to a text file with one gene per line. "
         "Mutually exclusive with --feature-types. If not specified, uses the default gene list for the selected --gb-type."
)
gene_filter_options.add_argument(
    "--feature-types",
    type=str,
    default=None,
    help="Comma-separated list of GenBank feature types to extract (e.g., 'CDS,rRNA,tRNA'). "
         "Mutually exclusive with --genes. Only applies to GenBank files. Default: 'CDS'."
)

pipeline_group = parser.add_argument_group(
    title="Pipeline Options",
    description="Options for controlling the analysis pipeline."
)
pipeline_group.add_argument(
    "--align",
    action="store_true",
    help="Perform multiple sequence alignment using MAFFT."
)
pipeline_group.add_argument(
    "--trimal",
    action="store_true",
    help="Perform alignment trimming using TrimAl. Requires --align."
)
pipeline_group.add_argument(
    "--iqtree",
    action="store_true",
    help="Perform phylogenetic analysis using IQ-TREE. Requires --trimal."
)
pipeline_group.add_argument(
    "--benchmark",
    action="store_true",
    help="Enable benchmarking of execution time."
)
pipeline_group.add_argument(
    "--allow-missing",
    action="store_true",
    help="Allow missing data in the supermatrix. Missing genes are filled with '?' characters. "
         "Without this flag, genes absent from any taxon are removed from all."
)
pipeline_group.add_argument(
    "--config",
    type=str,
    default=None,
    help="Path to a YAML configuration file for MAFFT, TrimAl, and IQ-TREE parameters. "
         "See tools_config.yaml for the default template."
)
pipeline_group.add_argument(
    "--overwrite",
    action="store_true",
    help="Overwrite existing output directories if they already exist."
)

if __name__ == "__main__":
    print(f"\n{'#'*70}\n")
    print(f"Version:     {__version__}")
    print(f"Status:      {__status__}")
    print(f"Authors:     {__authors__}")
    print(f"Maintainers: {__maintainers__}")
    print(f"Contact:     {__email__}")
    print(f"License:     {__license__}")
    print(f"GitHub:      {__github__}")
    print(f"\n{'#'*70}\n")
    
    tool_env = "splace_env"
    conda_env = os.getenv("CONDA_DEFAULT_ENV")
    in_container = os.path.exists("/.singularity.d") or os.getenv("SINGULARITY_CONTAINER") is not None

    if not in_container:
        if not shutil.which("conda"):
            logging.error("Conda is not installed or not found in PATH. Please install Conda and try again.")
            sys.exit(1)

        try:
            conda_envs = subprocess.run(
                ["conda", "env", "list"],
                capture_output=True,
                text=True,
                check=True
            ).stdout.splitlines()
        except subprocess.CalledProcessError:
            logging.error("Failed to list Conda environments.")
            sys.exit(1)

        if not any(tool_env in env for env in conda_envs):
            logging.error(f"The '{tool_env}' conda environment does not exist. Please create it using the provided environment.yml file.")
            sys.exit(1)

        if conda_env != tool_env:
            logging.warning(f"Activate the '{tool_env}' conda environment before running the {__tool__}. Current environment: {conda_env}. Run 'conda activate {tool_env}' and try again.")
            sys.exit(1)
    else:
        logging.info("Running inside a Singularity container.")

    args = parser.parse_args()

    if not args.input_dir and not args.ncbi_search_term:
        logging.error("Provide either --input_dir or --ncbi-search-term.")
        sys.exit(1)

    if args.ncbi_search_term and not args.ncbi_download_dir:
        logging.error("The argument --ncbi-search-term requires --ncbi-download-dir.")
        sys.exit(1)

    if args.ncbi_download_dir and not args.ncbi_search_term:
        logging.error("The argument --ncbi-download-dir requires --ncbi-search-term.")
        sys.exit(1)

    if args.download_only and not args.ncbi_search_term:
        logging.error("The argument --download-only requires --ncbi-search-term.")
        sys.exit(1)

    if args.ncbi_search_term and not (args.ncbi_complete or args.ncbi_partial):
        logging.error("The argument --ncbi-search-term requires at least one of --ncbi-complete or --ncbi-partial.")
        sys.exit(1)

    fasta_header_template = DEFAULT_FASTA_HEADER_CONFIG["template"]
    fasta_header_missing_value = DEFAULT_FASTA_HEADER_CONFIG["missing_value"]
    if not args.download_only:
        try:
            fasta_header_template, fasta_header_missing_value = load_fasta_header_config(args.fasta_header_config)
        except ValueError as exc:
            logging.error(str(exc))
            sys.exit(1)

        logging.info(f"Using FASTA header template from: {args.fasta_header_config}")

    # Parse gene filter: comma-separated string or text file
    genes_filter = None
    if args.genes:
        if os.path.isfile(args.genes):
            with open(args.genes, "r") as f:
                genes_filter = [line.strip() for line in f if line.strip()]
            logging.info(f"Loaded {len(genes_filter)} genes from file: {args.genes}")
        else:
            genes_filter = [g.strip() for g in args.genes.split(",") if g.strip()]
            logging.info(f"Using gene filter: {genes_filter}")

    # Parse feature types
    feature_types = None
    if args.feature_types:
        feature_types = [ft.strip() for ft in args.feature_types.split(",") if ft.strip()]
        logging.info(f"Using feature types: {feature_types}")

    # Load tool configuration
    tool_config = {
        "mafft": {"params": "--auto", "preserve_case": True, "timeout": 3600},
        "trimal": {"params": "-automated1", "timeout": 3600},
        "iqtree": {"bootstrap": 1000, "model": "MFP", "extra_args": ""},
    }
    if args.config:
        if not os.path.isfile(args.config):
            logging.error(f"Config file '{args.config}' not found.")
            sys.exit(1)
        with open(args.config, "r") as cfg_file:
            user_config = yaml.safe_load(cfg_file) or {}
        for section in tool_config:
            if section in user_config and isinstance(user_config[section], dict):
                tool_config[section].update(user_config[section])
        logging.info(f"Loaded tool configuration from: {args.config}")

    # Initialize Benchmark
    from splace.utils import Benchmark
    benchmark = Benchmark(output_path=os.path.join(args.output_dir, "benchmark.tsv"), enabled=args.benchmark)
        
    os.makedirs(name=args.output_dir, exist_ok=True, mode=0o755)

    input_directories = []
    if args.input_dir:
        input_directories.append(args.input_dir)

    downloaded_genbank_files = []
    if args.ncbi_search_term:
        from splace.ncbi_download import download_genbank_records_from_search

        benchmark.start("NCBI Download")
        try:
            downloaded_genbank_files = download_genbank_records_from_search(
                taxon_name=args.ncbi_search_term,
                output_dir=args.ncbi_download_dir,
                data_type=args.gb_type,
                include_complete=args.ncbi_complete,
                include_partial=args.ncbi_partial,
                refseq_only=args.ncbi_refseq_only,
                env_path=args.apis_env,
            )
        finally:
            benchmark.stop("NCBI Download")

        input_directories.append(args.ncbi_download_dir)

    if args.download_only:
        benchmark.save(
            input_dir=args.ncbi_download_dir or args.input_dir or os.getcwd(),
            file_count=len(downloaded_genbank_files),
        )
        logging.info(
            "NCBI retrieval stage completed. Add any outgroup files to the download directory before re-running SPLACE with extraction, alignment, trimming, or phylogeny flags."
        )
        sys.exit(0)

    try:
        genbank_files, fasta_files = collect_input_files(input_directories)
    except FileNotFoundError as exc:
        logging.error(str(exc))
        sys.exit(1)

    if downloaded_genbank_files:
        logging.info(
            f"Downloaded {len(downloaded_genbank_files)} GenBank file(s) from NCBI into '{args.ncbi_download_dir}'."
        )
    
    scanned_sources = ", ".join(input_directories)
    logging.info(msg=f"Found {len(fasta_files)} FASTA files and {len(genbank_files)} GenBank files in: {scanned_sources}")
    
    if genbank_files or fasta_files:
        from splace import extract_multiple_genbanks
        
        if args.trimal and not args.align:
            logging.error("The argument --trimal requires --align to be specified.")
            sys.exit(1)
        
        if args.iqtree and not args.trimal:
            logging.error("The argument --iqtree requires --trimal to be specified.")
            sys.exit(1)

        subdirs_to_check = ["markers_fasta", "aligned_markers", "phylogeny"]
        for subdir in subdirs_to_check:
            subdir_path = os.path.join(args.output_dir, subdir)
            if os.path.exists(subdir_path):
                if args.overwrite:
                    shutil.rmtree(subdir_path)
                    logging.info(f"Removed existing directory: {subdir_path}")
                else:
                    logging.error(f"Output directory '{subdir_path}' already exists. Use --overwrite to replace it.")
                    sys.exit(1)


        logging.info(msg=f"Found {len(genbank_files)} GenBank files and {len(fasta_files)} FASTA files. Processing...")
        try:
            converted_files = set()
            
            if genbank_files:
                benchmark.start("GenBank Extraction")
                gb_results = asyncio.run(
                    extract_multiple_genbanks(
                        genbank_files=genbank_files,
                        output_dir=os.path.join(args.output_dir, "markers_fasta"),
                        report_dir=args.output_dir,
                        data_type=args.gb_type,
                        max_concurrent=args.threads,
                        genes_filter=genes_filter,
                        feature_types=feature_types,
                        gbif_enabled=args.gbif,
                        header_template=fasta_header_template,
                        missing_value=fasta_header_missing_value,
                    )
                )
                benchmark.stop("GenBank Extraction")
                if gb_results:
                    converted_files.update(gb_results)

            if fasta_files:
                from splace import extract_multiple_fastas
                benchmark.start("FASTA Extraction")
                fasta_results = asyncio.run(
                    extract_multiple_fastas(
                        fasta_files=fasta_files,
                        output_dir=os.path.join(args.output_dir, "markers_fasta"),
                        data_type=args.gb_type,
                        max_concurrent=args.threads,
                        genes_filter=genes_filter
                    )
                )
                benchmark.stop("FASTA Extraction")
                if fasta_results:
                    converted_files.update(fasta_results)

            # Convert back to list for next steps
            converted_files = list(converted_files)

            if converted_files and args.align:
                from splace import align_multiple_files

                logging.info(msg=f"Starting alignment for {len(converted_files)} files...")
                benchmark.start("Alignment")
                aligned_files = asyncio.run(
                    align_multiple_files(
                        fasta_files=converted_files,
                        output_dir=os.path.join(args.output_dir, "aligned_markers"),
                        threads=args.threads,
                        mafft_params=tool_config["mafft"]["params"],
                        preserve_case=tool_config["mafft"]["preserve_case"],
                        timeout=tool_config["mafft"]["timeout"]
                    )
                )
                benchmark.stop("Alignment")

                if aligned_files and args.trimal:
                    from splace import trim_multiple_files

                    logging.info(msg=f"Starting trimming for {len(aligned_files)} files...")
                    benchmark.start("Trimming")
                    trimmed_files = asyncio.run(
                        trim_multiple_files(
                            alignment_files=aligned_files,
                            output_dir=os.path.join(args.output_dir, "trimmed_markers"),
                            trimal_params=tool_config["trimal"]["params"],
                            max_concurrent=args.threads,
                            timeout=tool_config["trimal"]["timeout"]
                        )
                    )
                    benchmark.stop("Trimming")

                    if trimmed_files and args.iqtree:
                        from splace import run_phylogeny_pipeline
                        
                        logging.info(msg=f"Starting phylogeny analysis...")
                        benchmark.start("Phylogeny")
                        run_phylogeny_pipeline(
                            trimmed_files=trimmed_files,
                            output_dir=os.path.join(args.output_dir, "phylogeny"),
                            threads=args.threads,
                            allow_missing=args.allow_missing,
                            bootstrap=tool_config["iqtree"]["bootstrap"],
                            model=tool_config["iqtree"]["model"],
                            extra_args=tool_config["iqtree"]["extra_args"],
                            report_dir=args.output_dir
                        )
                        benchmark.stop("Phylogeny")
            
            # Save Benchmark Statistics
            total_input_files = len(genbank_files) + len(fasta_files)
            benchmark.save(input_dir=args.input_dir, file_count=total_input_files)

        except Exception as e:
            logging.error(msg=f"An error occurred during processing: {e}")
    else:
        logging.error("No GenBank or FASTA files were found in the supplied input sources.")
        sys.exit(1)
    