#!/usr/bin/env python3
"""
download_fantom5_cage.py  -  start sites from FANTOM5 CAGE (hg38 v8 reprocessing)

Downloads three FANTOM5 files and writes one parsed table:

  datasets/raw/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz        peaks (BED9; recommended set)
  datasets/raw/human_phase1and2_CAGE_Peak_name.txt.gz            peak names, e.g. p1@TP53 (p1 = strongest peak of the gene)
  datasets/raw/hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz   peak -> EntrezGene ID, symbol
  datasets/raw/fantom5_cage_peaks_hg38_refseq.tsv                one row per peak, chromosome as RefSeq accession,
                                                                 representative TSS (1-based), gene ranks, Entrez IDs
  datasets/raw/fantom5_download_manifest.json                    URLs, sizes, MD5, counts

BED9 columns: chrom, start, end, peak ID, score, strand, thickStart (representative TSS,
0-based), thickEnd (= thickStart + 1), RGB.

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\download_fantom5_cage.py
    python scripts\\download_fantom5_cage.py --force     (re-download)
Standard library only.
"""
import argparse
import datetime as dt
import gzip
import hashlib
import json
import re
import ssl
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

BASES = ["https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/",
         "https://dbarchive.biosciencedbc.jp/data/fantom5/datafiles/reprocessed/hg38_v8/extra/"]   # mirror
FILES = {
    "bed": ("CAGE_peaks/", "hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"),
    "names": ("CAGE_peaks_annotation/", "human_phase1and2_CAGE_Peak_name.txt.gz"),
    "annot": ("CAGE_peaks_annotation/", "hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz"),
}
ASSEMBLY_REPORT = "GCF_000001405.40_GRCh38.p14_assembly_report.txt"
OUT = "fantom5_cage_peaks_hg38_refseq.tsv"


def md5_of(path, block=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def ssl_context():
    try:                                   # Windows Python often lacks the CA bundle: use certifi's
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def download(url, dest):
    print(f"Downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "HPSAD/1.0"})
    with urllib.request.urlopen(req, timeout=300, context=ssl_context()) as r, open(tmp, "wb") as out:
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
    m = {}
    with open(report_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[6] in ("na", ""):
                continue
            for alias in (f[0], f[4], f[6], f[9]):
                if alias and alias != "na":
                    m[alias] = f[6]
    return m


def id_keys(pid):
    """a peak ID may be 'old;new' in one file and one part in another: index all forms."""
    return [pid] + [p for p in pid.split(";") if p]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    raw = Path(a.root) / "datasets" / "raw"
    if not (raw / ASSEMBLY_REPORT).exists():
        sys.exit(f"Cannot find {raw / ASSEMBLY_REPORT}; run from C:\\Promoter\\HPSAD_v7 or pass --root.")
    used = {}
    for key, (sub, name) in FILES.items():
        if not a.force and (raw / name).exists():
            print(f"{name} already present; use --force to re-download")
            used[key] = "(already present)"
            continue
        for base in BASES:
            try:
                download(base + sub + name, raw / name)
                used[key] = base + sub + name
                break
            except Exception as e:                      # try the mirror next
                print(f"  failed: {type(e).__name__}: {e}")
        else:
            sys.exit(f"Could not download {name} from any source (see errors above).")
    chrom_map = load_chrom_map(raw / ASSEMBLY_REPORT)

    # peak names: p<rank>@<symbol>, comma-separated when a peak serves several genes
    ranks = defaultdict(dict)
    with gzip.open(raw / FILES["names"][1], "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            for nm in f[1].split(","):
                m = re.match(r"p(\d+)@(.+)", nm.strip())
                if m:
                    for k in id_keys(f[0]):
                        ranks[k][m.group(2)] = int(m.group(1))
    # Entrez IDs and symbols
    entrez = defaultdict(set)
    with gzip.open(raw / FILES["annot"][1], "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            for e in re.findall(r"(\d+)", f[3]):
                sym = f[7].strip()
                for k in id_keys(f[0]):
                    entrez[k].add((e, sym))

    n, unmapped, with_rank, with_entrez = 0, Counter(), 0, 0
    with gzip.open(raw / FILES["bed"][1], "rt", encoding="utf-8") as fh, \
            open(raw / OUT, "w", encoding="utf-8", newline="") as out:
        out.write("peak_id\tchrom\tucsc_chrom\tstrand\tpeak_start\tpeak_end\ttss\tscore\tgene_ranks\tentrez\n")
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            ref = chrom_map.get(f[0])
            if ref is None:
                unmapped[f[0]] += 1
                continue
            pid = f[3]
            rk, ez = {}, set()
            for k in id_keys(pid):
                rk.update(ranks.get(k, {}))
                ez |= entrez.get(k, set())
            with_rank += bool(rk)
            with_entrez += bool(ez)
            out.write("\t".join([pid, ref, f[0], f[5], str(int(f[1]) + 1), f[2], str(int(f[6]) + 1), f[4],
                                 ",".join(f"{s}:{r}" for s, r in sorted(rk.items(), key=lambda x: x[1])),
                                 ",".join(f"{e}:{s}" for e, s in sorted(ez))]) + "\n")
            n += 1
    summary = {"run_date": dt.datetime.now().isoformat(timespec="seconds"),
               "files": {k: {"url": used.get(k, ""), "bytes": (raw / name).stat().st_size, "md5": md5_of(raw / name)}
                         for k, (sub, name) in FILES.items()},
               "peaks_written": n, "peaks_with_gene_rank": with_rank, "peaks_with_entrez": with_entrez,
               "unmapped_chroms": dict(unmapped.most_common(20)), "parsed_tsv": OUT}
    json.dump(summary, open(raw / "fantom5_download_manifest.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))
    if with_rank == 0 or with_entrez == 0:
        print("WARNING: peak IDs did not join to the name or annotation file; paste the first lines of each file.")


if __name__ == "__main__":
    main()
