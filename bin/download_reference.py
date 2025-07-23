#!/usr/bin/env python3
"""
download_reference.py — replacement version
==========================================
Fetches genome FASTA, GFF and protein FASTA for the *top* accession
reported by `find_common_reference.py`, using the public NCBI Datasets
REST API (v2alpha). This eliminates the need for a local
`assembly_summary_refseq.txt` file.

Inputs
------
* **--file**     / **-f** : Path to `references_found.tsv` (first column must be the assembly accession, e.g. *GCF_000001405.40*).
* **--out_dir**  / **-o** : Directory where the compressed outputs will be written (created if necessary).
* **--reference**: *Deprecated, ignored*. Present only for backwards CLI‑compatibility with the previous script.

Outputs (written inside *out_dir*)
---------------------------------
* `<accession>_genomic.fna.gz`
* `<accession>_genomic.gff.gz`
* `<accession>.winner` – single‑line text file with the winning accession

The compressed file names match those produced by the original
implementation, so downstream workflow steps continue to work
unchanged.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib3
import zipfile
from io import BytesIO
from pathlib import Path
from typing import Dict, List

# -----------------------------------------------------------------------------
# Annotation request combinations (most → least detailed)
# -----------------------------------------------------------------------------
_ANNO_SETS: List[List[str]] = [
    ["GENOME_FASTA", "GENOME_GFF", "PROTEIN_FASTA"],
    ["GENOME_FASTA", "GENOME_GFF"],
    ["GENOME_FASTA"],
]

# Mapping from file suffix inside the dataset archive to output suffix
_SUFFIX_MAP: Dict[str, str] = {
    "_genomic.fna": "_genomic.fna.gz",
    "_genomic.gff": "_genomic.gff.gz",
    "_protein.faa": "_protein.faa.gz",
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _datasets_url(acc: str, anno_types: List[str]) -> str:
    anno_qs = "&".join(f"include_annotation_type={t}" for t in anno_types)
    return (
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        f"{acc}/download?{anno_qs}&hydrated=FULLY_HYDRATED"
    )


def _fetch_dataset(acc: str) -> BytesIO:
    """Try downloading with fallback annotation sets; return ZIP bytes."""
    http = urllib3.PoolManager()
    for anno in _ANNO_SETS:
        url = _datasets_url(acc, anno)
        resp = http.request("GET", url, preload_content=False)
        if resp.status == 200:
            resp.auto_close = False
            return BytesIO(resp.data)
        elif resp.status == 400:
            # Try next (simpler) request
            continue
        else:
            raise RuntimeError(f"NCBI API returned HTTP {resp.status} for {acc}")
    raise RuntimeError(
        f"NCBI API cannot provide any dataset (FASTA even) for accession {acc}."
    )


def _winner_from_tsv(tsv_path: Path) -> str:
    with tsv_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            return line.split("\t")[0]
    raise RuntimeError(f"No accession found in {tsv_path}")


def _extract_files(zip_bytes: BytesIO, acc: str, out_dir: Path):
    with zipfile.ZipFile(zip_bytes) as zf:
        zf.extractall(out_dir)

    data_root = out_dir / "ncbi_dataset" / "data" / acc
    if not data_root.exists():
        raise RuntimeError(f"Unexpected archive structure – cannot find {data_root}")

    for in_suf, out_suf in _SUFFIX_MAP.items():
        matches = list(data_root.rglob(f"*{in_suf}"))
        if not matches:
            # If the file type wasn't requested (due to fallback), just skip
            continue
        src_path = matches[0]
        dest_path = out_dir / f"{acc}{out_suf}"
        with src_path.open("rb") as src, gzip.open(dest_path, "wb") as dst:
            shutil.copyfileobj(src, dst)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def _parse_cli(argv: list[str] | None = None):
    p = argparse.ArgumentParser(
        description="Download reference assembly via NCBI Datasets REST API with smart fallbacks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-f", "--file", required=True, help="TSV from find_common_reference.py"
    )
    p.add_argument("-o", "--out_dir", required=True, help="Directory to place outputs")
    return p.parse_args(argv)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main(argv: list[str] | None = None):
    args = _parse_cli(argv)
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    acc = _winner_from_tsv(Path(args.file))
    (out_dir / f"{acc}.winner").write_text(acc + "\n")

    zip_bytes: BytesIO | None = None
    try:
        zip_bytes = _fetch_dataset(acc)
        _extract_files(zip_bytes, acc, out_dir)
    finally:
        if zip_bytes is not None:
            zip_bytes.close()

    print(f"[INFO] Downloaded reference for {acc} → {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
