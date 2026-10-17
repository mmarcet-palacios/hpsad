#!/usr/bin/env python3
"""
compare_tss_anchors.py  -  HPSAD fresh start, step 1: how good was the old anchor?

For every EPDnew promoter, finds the matching annotated gene (RefSeq GFF3, gene
and pseudogene features, primary chromosomes) and measures, in transcriptional
orientation:

  offset_gene_edge : EPDnew TSS minus the gene's 5' edge (the old HPSAD anchor)
  offset_nearest_tx: EPDnew TSS minus the nearest annotated transcript 5' end

Positive offset = EPDnew TSS lies downstream of the annotated site (inside the
gene); negative = upstream. 0 = exact agreement.

Outputs (datasets/qc/):
  tss_anchor_comparison.csv
  tss_anchor_comparison_summary.json
  tss_anchor_comparison.png   (only if matplotlib is installed)

Usage (from the project root):
    python scripts\\compare_tss_anchors.py
"""

import argparse
import bisect
import csv
import gzip
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path
from urllib.parse import unquote

GFF = "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
EPD = "Hs_EPDnew_006_hg38_refseq.tsv"
GENE_TYPES = {"gene", "pseudogene"}
BINS = [0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10**9]


def attrs(s):
    d = {}
    for part in s.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = unquote(v)
    return d


def five_prime(start, end, strand):
    return start if strand == "+" else end


def oriented(tss, site, strand):
    """distance from annotated site to EPD TSS, positive = downstream."""
    return (tss - site) if strand == "+" else (site - tss)


def load_gff(path):
    genes = {}                       # gene_id -> record
    by_symbol = defaultdict(list)
    tx_starts = defaultdict(set)     # gene_id -> transcript 5' ends
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[0].startswith("NC_0000"):
                continue
            chrom, ftype, start, end, strand = f[0], f[2], int(f[3]), int(f[4]), f[6]
            a = attrs(f[8])
            if ftype in GENE_TYPES:
                gid = a.get("ID")
                rec = {"id": gid, "chrom": chrom, "start": start, "end": end,
                       "strand": strand, "symbol": a.get("Name", a.get("gene", "")),
                       "biotype": a.get("gene_biotype", ftype),
                       "five": five_prime(start, end, strand)}
                genes[gid] = rec
                by_symbol[rec["symbol"]].append(rec)
                for syn in a.get("gene_synonym", "").split(","):
                    if syn:
                        by_symbol.setdefault("SYN:" + syn, []).append(rec)
            else:
                parent = a.get("Parent", "")
                if parent.startswith("gene-") and ftype not in ("exon", "CDS"):
                    tx_starts[parent].add(five_prime(start, end, strand))
    return genes, by_symbol, tx_starts


def match_gene(row, by_symbol, plus_index):
    chrom, strand, tss, sym = row["refseq_chrom"], row["strand"], int(row["tss_1based"]), row["gene_symbol"]
    for key, how in ((sym, "symbol"), ("SYN:" + sym, "synonym")):
        cands = [g for g in by_symbol.get(key, []) if g["chrom"] == chrom and g["strand"] == strand]
        if cands:
            return min(cands, key=lambda g: abs(oriented(tss, g["five"], strand))), how
    # fallback: nearest gene 5' end on same chromosome and strand within 5 kb
    fives, recs = plus_index.get((chrom, strand), ([], []))
    i = bisect.bisect_left(fives, tss)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(fives) and abs(fives[j] - tss) <= 5000:
            if best is None or abs(fives[j] - tss) < abs(best["five"] - tss):
                best = recs[j]
    return (best, "nearest_5kb") if best else (None, "unmatched")


def bin_label(x):
    ax = abs(x)
    for lo, hi in zip(BINS, BINS[1:]):
        if lo <= ax < hi:
            return "0" if hi == 1 else f"{lo}-{hi - 1}" if hi < 10**9 else f">={lo}"
    return "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    raw, qc = Path(a.root) / "datasets" / "raw", Path(a.root) / "datasets" / "qc"
    qc.mkdir(parents=True, exist_ok=True)

    print("Reading GFF3 ...")
    genes, by_symbol, tx_starts = load_gff(raw / GFF)
    idx = defaultdict(list)
    for g in genes.values():
        idx[(g["chrom"], g["strand"])].append((g["five"], g))
    plus_index = {k: ([x[0] for x in sorted(v, key=lambda t: t[0])],
                      [x[1] for x in sorted(v, key=lambda t: t[0])]) for k, v in idx.items()}
    print(f"  {len(genes):,} gene/pseudogene features on primary chromosomes")

    out_rows = []
    with open(raw / EPD, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g, how = match_gene(row, by_symbol, plus_index)
            tss, strand = int(row["tss_1based"]), row["strand"]
            r = {"epd_id": row["epd_id"], "gene_symbol": row["gene_symbol"],
                 "is_main_promoter": row["epd_id"].endswith("_1"),
                 "chrom": row["refseq_chrom"], "strand": strand, "epd_tss": tss,
                 "match": how, "gff_gene_id": "", "gff_biotype": "",
                 "gene_edge_5p": "", "offset_gene_edge": "", "offset_nearest_tx": ""}
            if g:
                txs = tx_starts.get(g["id"]) or {g["five"]}
                near = min(txs, key=lambda s: abs(tss - s))
                r.update({"gff_gene_id": g["id"], "gff_biotype": g["biotype"],
                          "gene_edge_5p": g["five"],
                          "offset_gene_edge": oriented(tss, g["five"], strand),
                          "offset_nearest_tx": oriented(tss, near, strand)})
            out_rows.append(r)

    with open(qc / "tss_anchor_comparison.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)

    def describe(rows, key):
        vals = sorted(r[key] for r in rows if r[key] != "")
        if not vals:
            return {}
        absv = sorted(abs(v) for v in vals)
        n = len(vals)
        pct = lambda p: absv[min(n - 1, int(p * n))]
        return {"n": n,
                "exact_0": sum(1 for v in vals if v == 0),
                "abs_median": pct(0.5), "abs_p75": pct(0.75), "abs_p90": pct(0.9),
                "epd_downstream_of_annotation": sum(1 for v in vals if v > 0),
                "epd_upstream_of_annotation": sum(1 for v in vals if v < 0),
                "abs_bins": dict(Counter(bin_label(v) for v in vals))}

    main_rows = [r for r in out_rows if r["is_main_promoter"]]
    summary = {
        "epd_promoters": len(out_rows),
        "match_method": dict(Counter(r["match"] for r in out_rows)),
        "matched_biotypes": dict(Counter(r["gff_biotype"] for r in out_rows if r["gff_biotype"])),
        "main_promoters_only": {"gene_edge": describe(main_rows, "offset_gene_edge"),
                                "nearest_transcript": describe(main_rows, "offset_nearest_tx")},
        "all_promoters": {"gene_edge": describe(out_rows, "offset_gene_edge"),
                          "nearest_transcript": describe(out_rows, "offset_nearest_tx")},
    }
    with open(qc / "tss_anchor_comparison_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
        for ax, key, title in ((axes[0], "offset_gene_edge", "vs gene 5' edge (old anchor)"),
                               (axes[1], "offset_nearest_tx", "vs nearest annotated transcript start")):
            v = [r[key] for r in main_rows if r[key] != "" and abs(r[key]) <= 500]
            ax.hist(v, bins=range(-500, 501, 10), color="#4C72B0")
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("EPDnew TSS offset (nt; + = downstream)")
        axes[0].set_ylabel("main EPDnew promoters")
        fig.tight_layout()
        fig.savefig(qc / "tss_anchor_comparison.png", dpi=200)
        print(f"Plot: {qc / 'tss_anchor_comparison.png'}")
    except ImportError:
        print("matplotlib not installed; skipped plot")


if __name__ == "__main__":
    main()
