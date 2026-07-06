import logging
import os
import time

from Bio import Entrez


DEFAULT_API_ENV_FILE = "apis.env"
DEFAULT_ENTREZ_EMAIL = "luan.rabelo@pq.itv.org"
FETCH_DELAY_DEFAULT = 0.5
FETCH_DELAY_WITH_KEY = 0.1


def load_api_settings(env_path=DEFAULT_API_ENV_FILE):
    """Load API settings from apis.env when available.

    Supported keys:
    - NCBI_API_KEY
    - NCBI_EMAIL
    """
    settings = {}
    if not env_path or not os.path.exists(env_path):
        return settings

    with open(env_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            settings[key.strip()] = value.strip().strip('"').strip("'")
    return settings


def find_api_settings(env_path=DEFAULT_API_ENV_FILE):
    """Resolve API settings from an explicit path, cwd, or environment variables."""
    candidate_paths = []
    if env_path:
        candidate_paths.append(env_path)
    candidate_paths.append(os.path.join(os.getcwd(), DEFAULT_API_ENV_FILE))

    settings = {}
    seen_paths = set()
    for candidate in candidate_paths:
        normalized = os.path.abspath(candidate)
        if normalized in seen_paths:
            continue
        seen_paths.add(normalized)
        if os.path.exists(normalized):
            settings.update(load_api_settings(normalized))
            settings["_source"] = normalized
            break

    if "NCBI_API_KEY" not in settings and os.getenv("NCBI_API_KEY"):
        settings["NCBI_API_KEY"] = os.getenv("NCBI_API_KEY", "").strip()
    if "NCBI_EMAIL" not in settings and os.getenv("NCBI_EMAIL"):
        settings["NCBI_EMAIL"] = os.getenv("NCBI_EMAIL", "").strip()

    return settings


def get_ncbi_credentials(api_key=None, email=None, env_path=DEFAULT_API_ENV_FILE):
    settings = find_api_settings(env_path=env_path)
    resolved_api_key = (api_key or settings.get("NCBI_API_KEY") or "").strip()
    resolved_email = (email or settings.get("NCBI_EMAIL") or DEFAULT_ENTREZ_EMAIL).strip()
    settings["NCBI_API_KEY"] = resolved_api_key
    settings["NCBI_EMAIL"] = resolved_email
    return settings


class EntrezRateLimiter:
    """Throttle Entrez requests to 10/s with an API key or 2/s otherwise."""

    def __init__(self, api_key=""):
        self.delay = FETCH_DELAY_WITH_KEY if api_key else FETCH_DELAY_DEFAULT
        self._last_request_at = 0.0

    def wait(self):
        if self._last_request_at:
            elapsed = time.monotonic() - self._last_request_at
            if elapsed < self.delay:
                time.sleep(self.delay - elapsed)
        self._last_request_at = time.monotonic()


def configure_entrez(api_key=None, email=None, env_path=DEFAULT_API_ENV_FILE):
    credentials = get_ncbi_credentials(api_key=api_key, email=email, env_path=env_path)
    Entrez.api_key = credentials.get("NCBI_API_KEY", "") or None
    Entrez.email = credentials.get("NCBI_EMAIL") or DEFAULT_ENTREZ_EMAIL
    return credentials


def create_entrez_session(api_key=None, email=None, env_path=DEFAULT_API_ENV_FILE):
    """Return configured credentials plus a shared rate limiter for one retrieval session."""
    credentials = configure_entrez(api_key=api_key, email=email, env_path=env_path)
    return credentials, EntrezRateLimiter(credentials.get("NCBI_API_KEY", ""))


def normalize_ncbi_taxon_name(taxon_name):
    normalized_taxon = str(taxon_name or "").strip()
    if not normalized_taxon:
        raise ValueError("NCBI taxon name cannot be empty.")
    return normalized_taxon


def get_ncbi_completeness_terms(include_complete=False, include_partial=False):
    completeness_terms = []
    if include_complete:
        completeness_terms.append("complete genome")
    if include_partial:
        completeness_terms.append("partial genome")

    if not completeness_terms:
        raise ValueError("At least one genome scope must be selected: include_complete or include_partial.")

    return completeness_terms


def get_ncbi_organelle_clause(data_type="mt"):
    if data_type == "cp":
        return "(chloroplast[Title] OR plastid[Title] OR chloroplast[Filter])"
    return "(mitochondrion[Title] OR mitochondrial[Title] OR mitochondrion[Filter])"


def build_ncbi_search_term(
    taxon_name,
    data_type="mt",
    include_complete=False,
    include_partial=False,
    refseq_only=False,
):
    """Build the final Entrez term used to search nucleotide genomes."""
    normalized_taxon = normalize_ncbi_taxon_name(taxon_name)
    completeness_terms = get_ncbi_completeness_terms(
        include_complete=include_complete,
        include_partial=include_partial,
    )

    if len(completeness_terms) == 1:
        completeness_clause = f'"{completeness_terms[0]}"[Title]'
    else:
        completeness_clause = "(" + " OR ".join(f'"{term}"[Title]' for term in completeness_terms) + ")"

    parts = [
        f'"{normalized_taxon}"[Organism]',
        completeness_clause,
        get_ncbi_organelle_clause(data_type=data_type),
        "biomol_genomic[PROP]",
    ]

    if refseq_only:
        parts.append("srcdb_refseq[PROP]")

    return " AND ".join(parts)


def search_ncbi_genome_ids(
    taxon_name,
    data_type="mt",
    include_complete=False,
    include_partial=False,
    refseq_only=False,
    retmax=200,
    api_key=None,
    email=None,
    env_path=DEFAULT_API_ENV_FILE,
    rate_limiter=None,
):
    """Search nucleotide IDs at NCBI using Bio.Entrez."""
    credentials = configure_entrez(api_key=api_key, email=email, env_path=env_path)
    if rate_limiter is None:
        rate_limiter = EntrezRateLimiter(credentials.get("NCBI_API_KEY", ""))
    rate_limiter.wait()

    with Entrez.esearch(
        db="nucleotide",
        term=build_ncbi_search_term(
            taxon_name=taxon_name,
            data_type=data_type,
            include_complete=include_complete,
            include_partial=include_partial,
            refseq_only=refseq_only,
        ),
        retmax=retmax,
        usehistory="n",
    ) as handle:
        payload = Entrez.read(handle)

    ids = payload.get("IdList", [])
    if len(ids) >= retmax:
        logging.warning(
            f"NCBI search reached the current limit of {retmax} records. Refine the search term if you need a smaller dataset."
        )
    return ids


def fetch_ncbi_summaries(ids, api_key=None, email=None, env_path=DEFAULT_API_ENV_FILE, rate_limiter=None):
    """Fetch esummary records for a list of NCBI nucleotide IDs."""
    if not ids:
        return {}

    credentials = configure_entrez(api_key=api_key, email=email, env_path=env_path)
    if rate_limiter is None:
        rate_limiter = EntrezRateLimiter(credentials.get("NCBI_API_KEY", ""))
    rate_limiter.wait()

    with Entrez.esummary(db="nucleotide", id=",".join(ids), retmode="xml") as handle:
        summary_records = Entrez.read(handle)

    return {str(record.get("Id", "")): record for record in summary_records}


def fetch_genbank_record(record_id, api_key=None, email=None, env_path=DEFAULT_API_ENV_FILE, rate_limiter=None):
    """Fetch one GenBank flatfile using Bio.Entrez."""
    credentials = configure_entrez(api_key=api_key, email=email, env_path=env_path)
    if rate_limiter is None:
        rate_limiter = EntrezRateLimiter(credentials.get("NCBI_API_KEY", ""))
    rate_limiter.wait()

    with Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text") as handle:
        text = handle.read()

    if "Error" in text and len(text) < 500:
        raise ValueError(f"NCBI returned an error for record {record_id}: {text.strip()}")
    return text


def _summary_accession(summary_record, fallback_id):
    return (
        summary_record.get("Caption")
        or summary_record.get("AccessionVersion")
        or summary_record.get("Id")
        or fallback_id
    )


def download_genbank_records_from_search(
    taxon_name,
    output_dir,
    data_type="mt",
    include_complete=False,
    include_partial=False,
    refseq_only=False,
    retmax=200,
    api_key=None,
    email=None,
    env_path=DEFAULT_API_ENV_FILE,
):
    """Search and download GenBank records from NCBI using Bio.Entrez.

    The function is intentionally small and reusable so it can be imported by
    future pipeline-oriented tooling.
    """
    if not output_dir:
        raise ValueError("An output directory is required to download GenBank files from NCBI.")

    os.makedirs(output_dir, exist_ok=True)
    credentials, rate_limiter = create_entrez_session(api_key=api_key, email=email, env_path=env_path)
    api_settings_source = credentials.get("_source")
    active_api_key = credentials.get("NCBI_API_KEY", "")
    active_email = credentials.get("NCBI_EMAIL", DEFAULT_ENTREZ_EMAIL)
    request_rate = 10 if active_api_key else 2

    if api_settings_source:
        logging.info(f"Loaded NCBI API settings from: {api_settings_source}")
    logging.info(
        f"Using Bio.Entrez for NCBI retrieval at up to {request_rate} request(s) per second with email '{active_email}'."
    )

    ids = search_ncbi_genome_ids(
        taxon_name=taxon_name,
        data_type=data_type,
        include_complete=include_complete,
        include_partial=include_partial,
        refseq_only=refseq_only,
        retmax=retmax,
        api_key=active_api_key,
        email=active_email,
        env_path=env_path,
        rate_limiter=rate_limiter,
    )

    if not ids:
        logging.warning("The NCBI search returned no records.")
        return []

    summaries = fetch_ncbi_summaries(
        ids,
        api_key=active_api_key,
        email=active_email,
        env_path=env_path,
        rate_limiter=rate_limiter,
    )
    downloaded_files = []

    for index, record_id in enumerate(ids, start=1):
        summary_record = summaries.get(str(record_id), {})
        accession = _summary_accession(summary_record, fallback_id=record_id)
        record_text = fetch_genbank_record(
            record_id,
            api_key=active_api_key,
            email=active_email,
            env_path=env_path,
            rate_limiter=rate_limiter,
        )
        output_path = os.path.join(output_dir, f"{accession}.gbk")
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write(record_text)

        downloaded_files.append(output_path)
        logging.info(f"Downloaded NCBI record {index}/{len(ids)}: {accession}")

    return downloaded_files