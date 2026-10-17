#!/usr/bin/env python3
"""
download_epdnew.py  -  HPSAD fresh start, step 0: start-site reference

Downloads the EPDnew human promoter collection (version 006, hg38) and writes:

  datasets/raw/Hs_EPDnew_006_hg38.sga.gz          original download (kept as-is)
  datasets/raw/Hs_EPDnew_006_hg38_refseq.tsv      parsed, one row per promoter,
                                                  chromosome mapped to NCBI accession
  datasets/raw/epdnew_download_manifest.json      URL, date, size, MD5, summary counts

EPDnew positions are experimentally supported transcription start sites
(1-based, single nucleotide). They replace gene-feature boundaries as the
position-zero anchor.

Usage (from C:\\Promoter\\HPSAD_GRCh38p14):
    python scripts\\download_epdnew.py
    python scripts\\download_epdnew.py --force     (re-download even if present)

Standard library only.
"""

import argparse
import datetime as dt
import gzip
import hashlib
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

SGA_URL = "https://epd.expasy.org/mga/hg38/epd/Hs_EPDnew_006_hg38.sga.gz"
SGA_NAME = "Hs_EPDnew_006_hg38.sga.gz"
ASSEMBLY_REPORT = "GCF_000001405.40_GRCh38.p14_assembly_report.txt"


def md5_of(path, block=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest):
    print(f"Downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "HPSAD/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as out:
        total = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
            total += len(chunk)
    tmp.replace(dest)
    print(f"  saved {dest.name} ({total:,} bytes)")


def load_chrom_map(report_path):
    """Map every naming style (UCSC chr1, GenBank CM..., RefSeq NC_...,
    sequence name '1') to the RefSeq accession, plus role and molecule type."""
    name_to_refseq, info = {}, {}
    with open(report_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            seq_name, role, mol, genbank, refseq, ucsc = f[0], f[1], f[3], f[4], f[6], f[9]
            if refseq in ("na", ""):
                continue
            info[refseq] = {"role": role, "molecule": mol, "ucsc": ucsc}
            for alias in (seq_name, genbank, refseq, ucsc):
                if alias and alias != "na":
                    name_to_refseq[alias] = refseq
    return name_to_refseq, info


def parse_sga(gz_path, name_to_refseq, info):
    """SGA columns: chrom, feature, position, strand, count, name[, ...]"""
    rows, unmapped = [], Counter()
    with gzip.open(gz_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                f = line.split()
            chrom, feature, pos, strand, count, name = f[:6]
            refseq = name_to_refseq.get(chrom)
            if refseq is None:
                unmapped[chrom] += 1
                continue
            rows.append({
                "epd_id": name,
                "gene_symbol": name.rsplit("_", 1)[0],
                "refseq_chrom": refseq,
                "ucsc_chrom": info[refseq]["ucsc"],
                "chrom_role": info[refseq]["role"],
                "tss_1based": int(pos),
                "strand": strand,
                "feature": feature,
                "count": count,
            })
    return rows, unmapped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="project root (HPSAD_GRCh38p14)")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    raw = Path(a.root) / "datasets" / "raw"
    if not (raw / ASSEMBLY_REPORT).exists():
        sys.exit(f"Cannot find {raw / ASSEMBLY_REPORT}. Run from the project root or pass --root.")

    sga = raw / SGA_NAME
    if a.force or not sga.exists():
        download(SGA_URL, sga)
    else:
        print(f"{sga.name} already present; use --force to re-download")

    name_to_refseq, info = load_chrom_map(raw / ASSEMBLY_REPORT)
    rows, unmapped = parse_sga(sga, name_to_refseq, info)

    out_tsv = raw / "Hs_EPDnew_006_hg38_refseq.tsv"
    cols = ["epd_id", "gene_symbol", "refseq_chrom", "ucsc_chrom", "chrom_role",
            "tss_1based", "strand", "feature", "count"]
    with open(out_tsv, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    ids = [r["epd_id"] for r in rows]
    genes = Counter(r["gene_symbol"] for r in rows)
    summary = {
        "source_url": SGA_URL,
        "downloaded_file": sga.name,
        "bytes": sga.stat().st_size,
        "md5": md5_of(sga),
        "run_date": dt.datetime.now().isoformat(timespec="seconds"),
        "promoters_parsed": len(rows),
        "unique_epd_ids": len(set(ids)),
        "duplicate_epd_ids": len(ids) - len(set(ids)),
        "unique_genes": len(genes),
        "genes_with_multiple_promoters": sum(1 for v in genes.values() if v > 1),
        "by_chrom_role": dict(Counter(r["chrom_role"] for r in rows)),
        "by_strand": dict(Counter(r["strand"] for r in rows)),
        "unmapped_chrom_names": dict(unmapped),
        "parsed_tsv": out_tsv.name,
    }
    with open(raw / "epdnew_download_manifest.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps({k: v for k, v in summary.items() if k != "source_url"}, indent=2))
    if unmapped:
        print("WARNING: some chromosome names did not map to RefSeq accessions (see manifest).")


if __name__ == "__main__":
    main()
