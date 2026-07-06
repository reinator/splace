import asyncio
import csv
import logging
import os
import re
import tempfile
import threading
import warnings

import requests

from Bio import BiopythonParserWarning, SeqIO
from SynGenes import SynGenes

from concurrent.futures import ThreadPoolExecutor

# Suppress BiopythonParserWarning globally — warnings.catch_warnings() is not
# thread-safe, so per-thread suppression inside ThreadPoolExecutor is unreliable.
warnings.filterwarnings("ignore", category=BiopythonParserWarning)

GENE_NAME_CACHE = {}
GBIF_CACHE = {}
FASTA_HEADER_FIELDS = {
    "accession",
    "authorship",
    "class",
    "family",
    "file_name",
    "genus",
    "kingdom",
    "order",
    "organism",
    "phylum",
    "species",
    "uid",
}
REPORT_METADATA_FIELDS = [
    "file_name",
    "accession",
    "organism",
    "gbif_status",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "authorship",
    "uid",
]
INVALID_SPECIES_FIELDS = ["file_name", "accession", "organism", "reason"]
MISSING_HEADER_FIELDS = ["file_name", "marker", "accession", "organism", "missing_fields", "header"]
HEADER_PLACEHOLDER_RE = re.compile(r"\{([a-zA-Z_][a-zA-Z0-9_]*)\}")
VALID_BINOMIAL_RE = re.compile(r"^[A-Z][A-Za-z.-]+\s+[a-z][A-Za-z.-]+$")

_uid_lock = threading.Lock()
_uid_counter = 0
_gbif_lock = threading.Lock()


def _next_uid():
    """Return the next unique 5-digit ID string, thread-safe."""
    global _uid_counter
    with _uid_lock:
        _uid_counter += 1
        return f"{_uid_counter:05d}"


sg = SynGenes(verbose=False)
# Workaround: SynGenes 1.0.6 builds log path without separator and may
# write to a read-only filesystem inside containers. Redirect to temp dir.
SynGenes.cwd_path = tempfile.gettempdir() + os.sep


def _extract_accession(record):
    accessions = getattr(record, "annotations", {}).get("accessions", []) or []
    if accessions:
        return accessions[0]
    return getattr(record, "id", "N/A")


def _split_organism_name(organism_name):
    parts = [part for part in re.split(r"\s+", str(organism_name or "").strip()) if part]
    genus = parts[0] if parts else ""
    species = parts[1] if len(parts) > 1 else ""
    return genus, species


def _extract_family_from_lineage(record):
    lineage = getattr(record, "annotations", {}).get("taxonomy", []) or []
    for entry in lineage:
        if re.search(r"(idae|aceae)$", str(entry), flags=re.IGNORECASE):
            return str(entry)
    return ""


def _is_valid_binomial_name(organism_name):
    return bool(VALID_BINOMIAL_RE.fullmatch(str(organism_name or "").strip()))


def _fetch_taxonomy_from_gbif(organism_name):
    cached = GBIF_CACHE.get(organism_name)
    if cached is not None:
        return cached

    try:
        response = requests.get(
            "https://api.gbif.org/v1/species",
            params={"name": organism_name},
            headers={"Accept": "application/json"},
            timeout=15,
        )
        response.raise_for_status()
        payload = response.json() or {}
        results = payload.get("results") or []
        if results:
            accepted = next((item for item in results if item.get("taxonomicStatus") == "ACCEPTED"), results[0])
            match = {
                "kingdom": accepted.get("kingdom", ""),
                "phylum": accepted.get("phylum", ""),
                "class": accepted.get("class", ""),
                "order": accepted.get("order", ""),
                "family": accepted.get("family", ""),
                "authorship": accepted.get("authorship", ""),
            }
        else:
            match = None
    except requests.RequestException as exc:
        logging.warning(f"GBIF lookup failed for '{organism_name}': {exc}")
        match = None

    with _gbif_lock:
        GBIF_CACHE[organism_name] = match
    return match


def _sanitize_header_value(value):
    text = str(value or "").strip()
    if not text:
        return ""
    return re.sub(r"\s+", "_", text)


def _build_fasta_header(template, metadata, missing_value):
    missing_fields = []

    def replace(match):
        field_name = match.group(1)
        value = _sanitize_header_value(metadata.get(field_name, ""))
        if value:
            return value
        missing_fields.append(field_name)
        return missing_value

    header = HEADER_PLACEHOLDER_RE.sub(replace, template)
    return header, list(dict.fromkeys(missing_fields))


def _extract_record_metadata(record, input_gb_file, record_uid, gbif_enabled):
    file_name = os.path.basename(input_gb_file)
    accession = _extract_accession(record)
    organism = getattr(record, "annotations", {}).get("organism", f"Unknown Species in {file_name}")
    genus, species = _split_organism_name(organism)
    metadata = {
        "accession": accession,
        "authorship": "",
        "class": "",
        "family": _extract_family_from_lineage(record),
        "file_name": file_name,
        "genus": genus,
        "kingdom": "",
        "order": "",
        "organism": organism,
        "phylum": "",
        "species": species,
        "uid": record_uid,
    }
    gbif_status = "not_requested"
    invalid_species = None

    if gbif_enabled:
        if _is_valid_binomial_name(organism):
            gbif_match = _fetch_taxonomy_from_gbif(organism)
            if gbif_match:
                metadata.update({key: value or metadata.get(key, "") for key, value in gbif_match.items()})
                gbif_status = "matched"
            else:
                gbif_status = "not_found"
        else:
            gbif_status = "invalid_name"
            invalid_species = {
                "file_name": file_name,
                "accession": accession,
                "organism": organism,
                "reason": "Expected a valid binomial in the format 'Genus species'.",
            }

    report_row = {field: metadata.get(field, "") for field in REPORT_METADATA_FIELDS}
    report_row["organism"] = organism
    report_row["gbif_status"] = gbif_status

    return metadata, report_row, invalid_species


def _write_tsv_report(file_path, fieldnames, rows):
    if not rows:
        return

    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_report_files(report_dir, metadata_rows, invalid_species_rows, missing_header_rows):
    metadata_path = os.path.join(report_dir, "metadata", "genbank_metadata.tsv")
    invalid_species_path = os.path.join(report_dir, "logs", "invalid_gbif_species.log")
    missing_header_path = os.path.join(report_dir, "logs", "fasta_header_missing_fields.log")

    _write_tsv_report(metadata_path, REPORT_METADATA_FIELDS, metadata_rows)
    _write_tsv_report(invalid_species_path, INVALID_SPECIES_FIELDS, invalid_species_rows)
    _write_tsv_report(missing_header_path, MISSING_HEADER_FIELDS, missing_header_rows)


async def genbank_to_fasta(**kwargs):
    """
    ## Convert GenBank file to FASTA format

    - Args:
        - `genbank_file` (**str**): Input GenBank file path
        - `output_path` (**str**): Output FASTA file path
        - `data_type` (**str**): Type of sequences to extract ("mt" for mitochondrial, "cp" for chloroplast)
        - `executor` (**ThreadPoolExecutor**): Executor for blocking I/O
        - `genes_filter` (**List[str]**, optional): List of specific gene names to extract. Overrides default genes_list.
        - `feature_types` (**List[str]**, optional): List of GenBank feature types to scan (e.g., ["CDS", "rRNA", "tRNA"]). Default: ["CDS"].

    - Returns:
        **dict**: Conversion outputs, metadata rows and report rows.
    """
    input_gb_file = kwargs.get("genbank_file", None)
    output_fasta_path = kwargs.get("output_path", None)
    data_type = kwargs.get("data_type", None)
    executor = kwargs.get("executor", None)
    genes_filter = kwargs.get("genes_filter", None)
    feature_types = kwargs.get("feature_types", None)
    gbif_enabled = kwargs.get("gbif_enabled", False)
    header_template = kwargs.get("header_template", "{organism}_{accession}_{uid}")
    missing_value = kwargs.get("missing_value", "?")

    if not input_gb_file or not os.path.exists(input_gb_file):
        logging.error(f"Input GenBank file ({input_gb_file}) not found.")
        return {
            "written_files": [],
            "metadata_rows": [],
            "invalid_species_rows": [],
            "missing_header_rows": [],
        }

    if output_fasta_path is None:
        output_fasta_path = os.getcwd()

    os.makedirs(output_fasta_path, exist_ok=True)

    def _process_genbank(input_gb_file, output_fasta_path, data_type, genes_filter, feature_types, gbif_enabled, header_template, missing_value):
        written_files = set()
        metadata_rows = []
        invalid_species_rows = []
        missing_header_rows = []

        if feature_types:
            target_feature_types = feature_types
        else:
            target_feature_types = ["CDS"]

        if genes_filter:
            genes_list = genes_filter
        elif data_type == "mt":
            genes_list = ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8"]
        elif data_type == "cp":
            genes_list = ["rbcL", "matK", "ndhF", "atpB", "psaA", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT"]
        else:
            genes_list = None

        best_per_gene = {}

        try:
            records = list(SeqIO.parse(input_gb_file, "genbank"))
            for record in records:
                record_uid = _next_uid()
                record_metadata, metadata_row, invalid_species = _extract_record_metadata(
                    record=record,
                    input_gb_file=input_gb_file,
                    record_uid=record_uid,
                    gbif_enabled=gbif_enabled,
                )
                metadata_rows.append(metadata_row)
                if invalid_species:
                    invalid_species_rows.append(invalid_species)

                header_body, header_missing_fields = _build_fasta_header(
                    template=header_template,
                    metadata=record_metadata,
                    missing_value=missing_value,
                )
                header = f">{header_body}"
                record_key = f"{record_metadata['organism']}_{record_metadata['accession']}_{record_uid}"

                for feature in record.features:
                    if feature.type not in target_feature_types:
                        continue

                    product = feature.qualifiers.get("product", [None])[0]
                    gene_qualifier = feature.qualifiers.get("gene", [None])[0]
                    raw_name = gene_qualifier or product

                    if not raw_name:
                        continue

                    if raw_name in GENE_NAME_CACHE:
                        gene_name = GENE_NAME_CACHE[raw_name]
                    else:
                        gene_name = sg.fix_gene_name(geneName=raw_name, type=str(data_type)) if data_type else raw_name
                        if gene_name and gene_name != "None":
                            GENE_NAME_CACHE[raw_name] = gene_name

                    if not gene_name or gene_name == "None":
                        logging.warning(
                            f"Could not standardize gene name for '{raw_name}' in file {os.path.basename(input_gb_file)}"
                        )
                        continue

                    if genes_list and gene_name not in genes_list:
                        continue

                    loc = feature.location
                    if hasattr(loc, "start") and hasattr(loc, "end") and int(loc.start) > int(loc.end):
                        logging.warning(
                            f"Skipping feature '{raw_name}' in {os.path.basename(input_gb_file)}: invalid location {loc} (start > end)"
                        )
                        continue

                    try:
                        seq_str = str(feature.extract(record.seq))
                    except Exception as exc:
                        logging.warning(
                            f"Skipping feature '{raw_name}' in {os.path.basename(input_gb_file)}: extraction failed ({exc})"
                        )
                        continue

                    seq_len = len(seq_str)
                    key = (gene_name, record_key)
                    if key not in best_per_gene or seq_len > best_per_gene[key]["length"]:
                        if key in best_per_gene:
                            logging.info(
                                f"Duplicate gene '{gene_name}' in {record_key}: keeping longer ({seq_len}bp over {best_per_gene[key]['length']}bp)"
                            )
                        best_per_gene[key] = {
                            "header": header,
                            "length": seq_len,
                            "marker": gene_name,
                            "missing_fields": header_missing_fields,
                            "organism": metadata_row["organism"],
                            "sequence": seq_str,
                            "accession": record_metadata["accession"],
                        }

            for (gene_name, _), entry in best_per_gene.items():
                output_file_path = os.path.join(output_fasta_path, f"{gene_name}.fasta")
                with open(output_file_path, "a+", encoding="utf-8") as output_file:
                    output_file.write(f"{entry['header']}\n{entry['sequence']}\n")
                written_files.add(output_file_path)

                if entry["missing_fields"]:
                    missing_header_rows.append(
                        {
                            "file_name": os.path.basename(input_gb_file),
                            "marker": gene_name,
                            "accession": entry["accession"],
                            "organism": entry["organism"],
                            "missing_fields": ",".join(entry["missing_fields"]),
                            "header": entry["header"][1:],
                        }
                    )

            if best_per_gene:
                logging.info(f"Converted {len(best_per_gene)} sequences from {os.path.basename(input_gb_file)} to FASTA")
            else:
                logging.warning(f"No matching sequences found in {os.path.basename(input_gb_file)}")

            return {
                "written_files": list(written_files),
                "metadata_rows": metadata_rows,
                "invalid_species_rows": invalid_species_rows,
                "missing_header_rows": missing_header_rows,
            }

        except Exception as exc:
            logging.error(f"Error converting {os.path.basename(input_gb_file)}: {exc}")
            return {
                "written_files": [],
                "metadata_rows": metadata_rows,
                "invalid_species_rows": invalid_species_rows,
                "missing_header_rows": missing_header_rows,
            }

    loop = asyncio.get_running_loop()
    if executor is None:
        with ThreadPoolExecutor() as local_executor:
            return await loop.run_in_executor(
                local_executor,
                _process_genbank,
                input_gb_file,
                output_fasta_path,
                data_type,
                genes_filter,
                feature_types,
                gbif_enabled,
                header_template,
                missing_value,
            )

    return await loop.run_in_executor(
        executor,
        _process_genbank,
        input_gb_file,
        output_fasta_path,
        data_type,
        genes_filter,
        feature_types,
        gbif_enabled,
        header_template,
        missing_value,
    )


async def extract_multiple_genbanks(**kwargs):
    """
    ## Convert multiple GenBank files to FASTA with concurrency control

    - Args:
        - `genbank_files` (**List[str]**): List of GenBank file paths
        - `output_dir` (**str**): Output directory for FASTA files
        - `data_type` (**str**): Type filter ("mt" or "cp")
        - `max_concurrent` (**int**): Maximum concurrent conversions (default: 5)
        - `genes_filter` (**List[str]**, optional): List of specific gene names to extract
        - `feature_types` (**List[str]**, optional): List of GenBank feature types to scan

    - Returns:
        **List[str]**: List of successfully converted file paths
    """
    list_gb_files = kwargs.get("genbank_files", [])
    output_directory = kwargs.get("output_dir", "markers_fasta")
    report_directory = kwargs.get("report_dir") or output_directory
    data_type = kwargs.get("data_type", None)
    max_concurrent = kwargs.get("max_concurrent", 5)
    genes_filter = kwargs.get("genes_filter", None)
    feature_types = kwargs.get("feature_types", None)
    gbif_enabled = kwargs.get("gbif_enabled", False)
    header_template = kwargs.get("header_template", "{organism}_{accession}_{uid}")
    missing_value = kwargs.get("missing_value", "?")

    if not list_gb_files:
        logging.warning("No GenBank files provided for conversion")
        return []

    if not os.path.exists(output_directory):
        logging.info(f"Creating output directory: {output_directory}")
        os.makedirs(output_directory, exist_ok=True, mode=0o755)

    semaphore = asyncio.Semaphore(max_concurrent)
    executor = ThreadPoolExecutor(max_workers=max_concurrent)

    async def convert_single_file(gb_file):
        async with semaphore:
            return await genbank_to_fasta(
                genbank_file=gb_file,
                output_path=output_directory,
                data_type=data_type,
                executor=executor,
                genes_filter=genes_filter,
                feature_types=feature_types,
                gbif_enabled=gbif_enabled,
                header_template=header_template,
                missing_value=missing_value,
            )

    tasks = [convert_single_file(gb_file) for gb_file in list_gb_files]
    logging.info(f"Starting conversion of {len(list_gb_files)} GenBank files")

    results = await asyncio.gather(*tasks, return_exceptions=True)
    executor.shutdown(wait=False)

    converted_files = set()
    metadata_rows = []
    invalid_species_rows = []
    missing_header_rows = []

    for index, result in enumerate(results):
        if isinstance(result, Exception):
            logging.error(f"Error converting {list_gb_files[index]}: {result}")
            continue

        if not isinstance(result, dict):
            continue

        converted_files.update(result.get("written_files", []))
        metadata_rows.extend(result.get("metadata_rows", []))
        invalid_species_rows.extend(result.get("invalid_species_rows", []))
        missing_header_rows.extend(result.get("missing_header_rows", []))

    _write_report_files(
        report_dir=report_directory,
        metadata_rows=metadata_rows,
        invalid_species_rows=invalid_species_rows,
        missing_header_rows=missing_header_rows,
    )

    logging.info(f"Conversion completed. {len(converted_files)} unique gene files generated/updated.")
    return list(converted_files)
