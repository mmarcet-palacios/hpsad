#!/usr/bin/env python3
"""
match_cage_anchors_v7.py  -  give every gene its strongest FANTOM5 CAGE start site

For each gene and pseudogene in the GFF, pick one CAGE peak as its start site:
  1. entrez    peaks associated with the gene's GeneID (FANTOM5 annotation); the gene's
               lowest-rank peak (p1 = strongest across FANTOM5 samples) wins, ties -> nearer the 5' edge
  2. symbol    if no GeneID match: peaks named p<rank>@<gene symbol>
  3. geometric if neither: an unassigned peak (no p<rank>@gene name) on the same strand within
               +/-500 nt of the annotated 5' edge, the nearest one
A chosen peak must lie on the gene's strand, between 1,000 nt upstream of the annotated 5'
edge and the gene's 3' end. Peaks assigned only to OTHER genes are never used.

Validation: for genes that also have an EPDnew main promoter, the offset CAGE - EPDnew.

Usage (from C:\\Promoter\\HPSAD_v7, after download_fantom5_cage.py and compare_tss_anchors.py):
    python scripts\\match_cage_anchors_v7.py
Outputs: datasets/qc/cage_anchor_match.csv, cage_anchor_match_summary.json
"""
import argparse
import bisect
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

GFF = "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
PEAKS = "fantom5_cage_peaks_hg38_refseq.tsv"
EPD_MATCH = "tss_anchor_comparison.csv"
UP, GEO = 1000, 500


def gene_class(ftype, biotype):
    if ftype == "pseudogene" or "pseudogene" in biotype:
        return "pseudogene"
    return "protein_coding" if biotype == "protein_coding" else "non_coding"


def load_genes(path):
    genes = {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[0].startswith("NC_0000") or f[2] not in ("gene", "pseudogene"):
                continue
            a = dict(p.split("=", 1) for p in f[8].split(";") if "=" in p)
            gid = a.get("ID", "")
            s, e, strand = int(f[3]), int(f[4]), f[6]
            bt = unquote(a.get("gene_biotype", f[2]))
            m = re.search(r"GeneID:(\d+)", unquote(a.get("Dbxref", "")))
            genes[gid] = dict(chrom=f[0], start=s, end=e, strand=strand, symbol=unquote(a.get("Name", "")),
                              biotype=bt, gclass=gene_class(f[2], bt), entrez=m.group(1) if m else "",
                              five=s if strand == "+" else e)
    return genes


def oriented(pos, five, strand):
    return pos - five if strand == "+" else five - pos


def bin_abs(x):
    x = abs(x)
    for lim, lab in ((0, "0"), (10, "1-10"), (50, "11-50"), (100, "51-100"), (500, "101-500"), (5000, "501-5000")):
        if x <= lim:
            return lab
    return ">5000"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    root = Path(a.root)
    raw, qc = root / "datasets" / "raw", root / "datasets" / "qc"
    print("Reading GFF3 ...")
    genes = load_genes(raw / GFF)
    peaks = defaultdict(list)
    with open(raw / PEAKS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ranks = dict((x.rsplit(":", 1)[0], int(x.rsplit(":", 1)[1])) for x in r["gene_ranks"].split(",") if x)
            ez = defaultdict(set)
            for x in r["entrez"].split(","):
                if x:
                    e, s = x.split(":", 1)
                    ez[e].add(s)
            peaks[(r["chrom"], r["strand"])].append((int(r["tss"]), r["peak_id"], ranks, dict(ez)))
    for k in peaks:
        peaks[k].sort(key=lambda x: x[0])
    pos_index = {k: [p[0] for p in v] for k, v in peaks.items()}
    print(f"  {len(genes):,} genes; {sum(len(v) for v in peaks.values()):,} CAGE peaks")

    epd = {}
    with open(qc / EPD_MATCH, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["is_main_promoter"] == "True" and r["gff_gene_id"] in genes:
                epd[r["gff_gene_id"]] = int(r["epd_tss"])

    rows = []
    for gid, g in genes.items():
        key = (g["chrom"], g["strand"])
        if key not in peaks:
            continue
        glen = g["end"] - g["start"] + 1
        lo, hi = (g["five"] - UP, g["five"] + glen) if g["strand"] == "+" else (g["five"] - glen, g["five"] + UP)
        i0, i1 = bisect.bisect_left(pos_index[key], lo), bisect.bisect_right(pos_index[key], hi)
        cand = peaks[key][i0:i1]
        best, method = None, ""
        for meth in ("entrez", "symbol", "geometric"):
            pool = []
            for tss, pid, ranks, ez in cand:
                off = oriented(tss, g["five"], g["strand"])
                if meth == "entrez" and g["entrez"] in ez:
                    r = min([ranks.get(s, 999) for s in ez[g["entrez"]]] + [ranks.get(g["symbol"], 999)])
                    pool.append((r, abs(off), tss, pid))
                elif meth == "symbol" and g["symbol"] in ranks:
                    pool.append((ranks[g["symbol"]], abs(off), tss, pid))
                elif meth == "geometric" and not ranks and not ez and abs(off) <= GEO:
                    pool.append((999, abs(off), tss, pid))
            if pool:
                best, method = min(pool), meth
                break
        if best is None:
            continue
        rank, _, tss, pid = best
        e = epd.get(gid)
        rows.append(dict(gene_id=gid, symbol=g["symbol"], gene_class=g["gclass"], biotype=g["biotype"],
                         chrom=g["chrom"], strand=g["strand"], gene_edge_5p=g["five"], cage_peak=pid,
                         cage_tss=tss, cage_rank=rank if rank < 999 else "", method=method,
                         offset_gene_edge=oriented(tss, g["five"], g["strand"]),
                         epd_main_tss=e or "", offset_cage_vs_epd=oriented(tss, e, g["strand"]) if e else ""))

    with open(qc / "cage_anchor_match.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    by_class = Counter(g["gclass"] for g in genes.values())
    matched = Counter(r["gene_class"] for r in rows)
    no_epd = [r for r in rows if not r["epd_main_tss"]]
    val = [int(r["offset_cage_vs_epd"]) for r in rows if r["offset_cage_vs_epd"] != ""]
    order = ["0", "1-10", "11-50", "51-100", "101-500", "501-5000", ">5000"]
    summary = {
        "genes_by_class": dict(by_class),
        "matched_by_class": dict(matched),
        "matched_by_method": dict(Counter(r["method"] for r in rows)),
        "matched_without_epdnew_by_class": dict(Counter(r["gene_class"] for r in no_epd)),
        "rank1_share_without_epdnew": round(sum(1 for r in no_epd if r["cage_rank"] == 1) / max(len(no_epd), 1), 3),
        "offset_cage_vs_gene_edge_without_epdnew": {c: dict((b, n) for b, n in sorted(
            Counter(bin_abs(r["offset_gene_edge"]) for r in no_epd if r["gene_class"] == c).items(),
            key=lambda x: order.index(x[0]))) for c in ("protein_coding", "non_coding", "pseudogene")},
        "validation_cage_vs_epdnew": {"genes": len(val),
                                      "abs_bins": dict((b, n) for b, n in sorted(Counter(bin_abs(v) for v in val).items(),
                                                                                key=lambda x: order.index(x[0]))),
                                      "within_10nt": round(sum(abs(v) <= 10 for v in val) / max(len(val), 1), 3),
                                      "within_50nt": round(sum(abs(v) <= 50 for v in val) / max(len(val), 1), 3)},
    }
    json.dump(summary, open(qc / "cage_anchor_match_summary.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
