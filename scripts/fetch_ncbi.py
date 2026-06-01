#!/usr/bin/env python3
"""
Fetch sequences and metadata from NCBI Entrez for use in the flexpipe pipeline.

Outputs are in Pathoplexus-compatible TSV format so that merge_local_sequences.py
and curate.py can process Pathoplexus and NCBI data uniformly.

Column mapping (GenBank → Pathoplexus):
    accession.version     → accessionVersion
    collection_date       → sampleCollectionDate
    country (parsed)      → geoLocCountry, geoLocAdmin1, geoLocAdmin2
    host                  → hostNameCommon
    authors               → authors
    fixed "OPEN"          → dataUseTerms
    fixed ""              → lineage   (Nextclade fills this later)
    fixed "NCBI"          → source
"""
import argparse
import os
import re
import socket
import sys
import time
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError

import pandas as pd
import yaml
from Bio import Entrez, SeqIO

# INSDC-compliant date patterns — fully anchored so ranges like
# 'Feb-2015/Apr-2015' are rejected rather than partially matched
_DATE_RE = re.compile(
    r"^("
    r"\d{4}(-\d{2}(-\d{2})?)?"   # YYYY, YYYY-MM, YYYY-MM-DD
    r"|\d{2}-[A-Za-z]{3}-\d{4}"  # DD-Mon-YYYY  (e.g. 15-Jan-2020)
    r"|[A-Za-z]{3}-\d{4}"        # Mon-YYYY     (e.g. Jan-2020)
    r")$"
)


def clean_date(raw: str) -> str:
    """Return raw if it looks like an INSDC date, otherwise empty string."""
    raw = raw.strip()
    return raw if _DATE_RE.match(raw) else ""


BATCH_SIZE = 200
DELAY_SEC  = 0.4   # ≤ 3 req/s without API key, ≤ 10 with API key


# ── parsing helpers ──────────────────────────────────────────────────────────

def parse_country_field(raw):
    """Parse NCBI country field 'Country: Division, Location' → (country, div, loc)."""
    country = division = location = ""
    if not raw:
        return country, division, location
    if ":" in raw:
        country, sub = raw.split(":", 1)
        country = country.strip()
        sub = sub.strip()
        if "," in sub:
            parts = [p.strip() for p in sub.split(",", 1)]
            division, location = parts[0], parts[1]
        else:
            division = sub
    else:
        country = raw.strip()
    return country, division, location


def parse_gb_record(rec):
    """Extract metadata dict from a BioPython GenBank SeqRecord."""
    host            = ""
    raw_country     = ""
    collection_date = ""
    authors         = ""

    for feature in rec.features:
        if feature.type == "source":
            host            = feature.qualifiers.get("host",            [""])[0]
            # INSDC migrated from "country" to "geo_loc_name" qualifier (~2023)
            raw_country     = (feature.qualifiers.get("geo_loc_name") or
                               feature.qualifiers.get("country") or [""])[0]
            collection_date = clean_date(
                feature.qualifiers.get("collection_date", [""])[0]
            )
            break

    for ref in rec.annotations.get("references", []):
        if ref.authors and not authors:
            first = ref.authors.split(",")[0].strip()
            authors = f"{first} et al"

    country, division, location = parse_country_field(raw_country)

    return {
        "accessionVersion":     rec.id,
        "sampleCollectionDate": collection_date,
        "geoLocCountry":        country,
        "geoLocAdmin1":         division,
        "geoLocAdmin2":         location,
        "hostNameCommon":       host,
        "authors":              authors,
        "dataUseTerms":         "OPEN",
        "lineage":              "",
        "source":               "NCBI",
    }


# ── NCBI search + fetch ───────────────────────────────────────────────────────

def search_ncbi(taxid, min_length, max_length, min_date=None, extra_term=None):
    """Return (count, webenv, query_key) for the NCBI result set."""
    query = f"txid{taxid}[Organism] {min_length}:{max_length}[SLEN]"
    if min_date:
        # NCBI PDAT requires YYYY/MM/DD format
        ncbi_date = str(min_date).replace("-", "/")
        query += f" {ncbi_date}:3000/12/31[PDAT]"
    if extra_term:
        query += f" {extra_term}"
    print(f"NCBI query: {query}", flush=True)

    handle = Entrez.esearch(
        db="nucleotide", term=query, idtype="acc", usehistory="y"
    )
    result = Entrez.read(handle)
    handle.close()

    count = int(result["Count"])
    print(f"Found {count} records on NCBI.", flush=True)
    return count, result["WebEnv"], result["QueryKey"]


def iter_records(count, webenv, query_key):
    """Yield BioPython SeqRecords from NCBI server history in batches, with retries."""
    _transient = (IncompleteRead, HTTPError, URLError, socket.error, OSError)
    for start in range(0, count, BATCH_SIZE):
        end = min(start + BATCH_SIZE, count)
        print(f"  Fetching records {start + 1}–{end} / {count} ...", flush=True)
        for attempt in range(1, 6):
            handle = None
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb", retmode="text",
                    retstart=start, retmax=BATCH_SIZE,
                    webenv=webenv, query_key=query_key,
                )
                records = list(SeqIO.parse(handle, "gb"))
                handle.close()
                break
            except _transient as exc:
                if handle is not None:
                    try:
                        handle.close()
                    except Exception:
                        pass
                if attempt < 5:
                    wait = 5 * attempt
                    print(f"    Network error ({exc.__class__.__name__}), retry {attempt}/5 in {wait}s...", flush=True)
                    time.sleep(wait)
                else:
                    raise
        for rec in records:
            yield rec
        time.sleep(DELAY_SEC)


# ── main ──────────────────────────────────────────────────────────────────────

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch sequences and metadata from NCBI for the flexpipe pipeline."
    )
    parser.add_argument("--config",           required=True,
                        help="Path to config/config.yaml")
    parser.add_argument("--metadata-output",  required=True,
                        help="Output TSV path (Pathoplexus-compatible format)")
    parser.add_argument("--sequences-output", required=True,
                        help="Output FASTA path")
    args = parser.parse_args()

    cfg  = load_config(args.config)
    ncbi = cfg.get("ncbi", {})
    sub  = cfg.get("subsampling", {})

    taxid       = ncbi.get("taxid")
    genome_size = ncbi.get("genome_size")
    email       = ncbi.get("email", "") or "pipeline@example.com"
    api_key     = ncbi.get("api_key", "") or None
    min_frac    = float(ncbi.get("min_length", 0.7))
    max_frac    = float(ncbi.get("max_length", 1.1))

    if not taxid:
        sys.exit("ERROR: ncbi.taxid is required in config.yaml")
    if not genome_size:
        sys.exit("ERROR: ncbi.genome_size is required in config.yaml")

    min_length = int(genome_size * min_frac)
    max_length = int(genome_size * max_frac)

    # min_date: explicit ncbi.min_date falls back to subsampling.min_year
    min_date = ncbi.get("min_date") or None
    if not min_date:
        min_year = sub.get("min_year")
        if min_year:
            min_date = str(min_year)

    extra_term = ncbi.get("extra_search_term") or None

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    os.makedirs(os.path.dirname(os.path.abspath(args.metadata_output)),  exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.sequences_output)), exist_ok=True)

    print(f"Length filter: {min_length}–{max_length} bp  |  min_date: {min_date or 'none'}")
    if extra_term:
        print(f"Extra search term: {extra_term}")

    count, webenv, query_key = search_ncbi(taxid, min_length, max_length, min_date, extra_term)

    # empty outputs so Snakemake never fails on a 0-result query
    if count == 0:
        pd.DataFrame(columns=[
            "accessionVersion", "sampleCollectionDate", "geoLocCountry",
            "geoLocAdmin1", "geoLocAdmin2", "hostNameCommon", "authors",
            "dataUseTerms", "lineage", "source",
        ]).to_csv(args.metadata_output, sep="\t", index=False)
        open(args.sequences_output, "w").close()
        print("No records found — empty outputs written.")
        return

    rows  = []
    n_seq = 0

    with open(args.sequences_output, "w") as fa:
        for rec in iter_records(count, webenv, query_key):
            rows.append(parse_gb_record(rec))
            fa.write(f">{rec.id}\n{rec.seq}\n")
            n_seq += 1

    df = pd.DataFrame(rows)

    # Self-healing: NCBI sometimes indexes new bulk submissions without geo_loc_name
    # populated yet. If >5 % of records have empty country, re-fetch those records
    # individually to pick up any qualifier that arrived after the initial batch fetch.
    empty_mask = df["geoLocCountry"] == ""
    n_empty = int(empty_mask.sum())
    if n_empty > 0 and n_empty / len(df) > 0.05:
        print(f"\n[repair] {n_empty}/{len(df)} records have empty geoLocCountry — re-fetching...",
              flush=True)
        empty_accs = df.loc[empty_mask, "accessionVersion"].tolist()
        _transient = (IncompleteRead, HTTPError, URLError, socket.error, OSError)
        for start in range(0, len(empty_accs), BATCH_SIZE):
            batch = empty_accs[start:start + BATCH_SIZE]
            for attempt in range(1, 6):
                handle = None
                try:
                    handle = Entrez.efetch(
                        db="nucleotide", id=",".join(batch),
                        rettype="gb", retmode="text",
                    )
                    fixed = {r.id: parse_gb_record(r)
                             for r in SeqIO.parse(handle, "gb")}
                    handle.close()
                    break
                except _transient as exc:
                    if handle:
                        try: handle.close()
                        except Exception: pass
                    if attempt < 5:
                        time.sleep(5 * attempt)
                    else:
                        fixed = {}
            for acc, meta in fixed.items():
                if meta["geoLocCountry"]:
                    idx = df.index[df["accessionVersion"] == acc]
                    for col in ("geoLocCountry", "geoLocAdmin1", "geoLocAdmin2"):
                        df.loc[idx, col] = meta[col]
            time.sleep(DELAY_SEC)
        filled = int((df["geoLocCountry"] != "").sum()) - (len(df) - n_empty)
        print(f"[repair] filled {filled} country values", flush=True)

    df.to_csv(args.metadata_output, sep="\t", index=False)

    print(f"\nMetadata:  {len(df)} records → {args.metadata_output}")
    print(f"Sequences: {n_seq} records  → {args.sequences_output}")


if __name__ == "__main__":
    main()
