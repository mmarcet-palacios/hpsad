#!/usr/bin/env python3
"""
build_cohorts_v7.py  -  HPSAD fresh start, step 2a: anchors, segments and cohorts
(coordinates only; sequences are added by extract_cohort_sequences_v7.py)

ANCHORS (position zero)
  * Genes matched to EPDnew: the main EPDnew promoter (_1) is the primary anchor;
    other EPDnew promoters of the same gene are alternative anchors.
  * With --cage: genes without EPDnew but with a FANTOM5 CAGE start site
    (match_cage_anchors_v7.py) take that site as the primary anchor; their
    annotated transcript 5' ends (gene edge included) become alternative anchors.
  * All other genes and pseudogenes: the gene 5' edge is the primary anchor;
    other distinct annotated transcript 5' ends are alternative anchors.
  anchor_source records which (EPDnew / CAGE / annotation).

SEGMENT RULES (upstream of the anchor, up to 2,000 nt)
  R1  The segment may not overlap any OTHER gene or pseudogene, on either strand.
      If one is in the way, the segment is cut where it meets that gene.
  R2  ALTERNATIVE anchors only: the segment may not contain intron sequence of
      its OWN gene (intron = inside the gene span but outside every annotated exon
      of that gene). Primary anchors ignore their own gene.
  C4 anchors (start site inside another gene): R1 uses the other genes' EXONS
      only, so host-gene intron sequence is allowed.
  N-containing segments are removed later, at sequence extraction.

COHORTS
  C1   primary anchor, full 2,000 nt
  C2   primary anchor, truncated (1-1,999 nt)
  C2b  as C2, and an opposite-strand gene starts within 1,000 nt upstream (head-to-head)
  C3   alternative anchor, segment >= MIN_ALT nt after R1 and R2
  C3x  alternative anchor failing (a): internal to its own gene (catalogued, kept for later)
  C4   anchor lies inside another gene or pseudogene (nested); segment cut only at
       other genes' exons. Analysed as its own cohort, never mixed with C1/C2
  C0   primary anchor with no usable upstream sequence (catalogued)
Every row carries gene_class (protein_coding / non_coding / pseudogene) and anchor_source.

Outputs (datasets/cohorts/):
  hpsad_anchors_v7.tsv            one row per anchor, all cohorts
  hpsad_anchors_v7_summary.json   counts by cohort x gene_class x anchor_source

Usage (from the project root):  python scripts\\build_cohorts_v7.py [--cage]
"""

import argparse
import bisect
import csv
import gzip
import json
from collections import defaultdict, Counter
from pathlib import Path
from urllib.parse import unquote

GFF = "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
EPD = "Hs_EPDnew_006_hg38_refseq.tsv"
EPD_MATCH = "tss_anchor_comparison.csv"          # from compare_tss_anchors.py
L_MAX = 2000
H2H_WINDOW = 1000
MIN_ALT = 50
GENE_TYPES = {"gene", "pseudogene"}
SKIP_TYPES = {"exon", "CDS"}


def attrs(s):
    d = {}
    for part in s.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = unquote(v)
    return d


def gene_class(ftype, biotype):
    if ftype == "pseudogene" or "pseudogene" in biotype:
        return "pseudogene"
    return "protein_coding" if biotype == "protein_coding" else "non_coding"


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


# ── GFF ─────────────────────────────────────────────────────────────────────
def load_gff(path):
    genes, parent_of, tx5, exons_raw = {}, {}, defaultdict(set), defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[0].startswith("NC_0000"):
                continue
            chrom, ftype, s, e, strand = f[0], f[2], int(f[3]), int(f[4]), f[6]
            a = attrs(f[8])
            fid, par = a.get("ID", ""), a.get("Parent", "").split(",")[0]
            if ftype in GENE_TYPES:
                bt = a.get("gene_biotype", ftype)
                genes[fid] = dict(id=fid, chrom=chrom, start=s, end=e, strand=strand,
                                  symbol=a.get("Name", a.get("gene", "")), biotype=bt,
                                  gclass=gene_class(ftype, bt),
                                  five=s if strand == "+" else e)
            elif ftype == "exon":
                exons_raw[par].append((s, e))
            elif ftype not in SKIP_TYPES and fid and par:
                parent_of[fid] = par
                if par.startswith("gene-"):
                    tx5[par].add(s if strand == "+" else e)

    def to_gene(x):
        seen = 0
        while x and not x.startswith("gene-") and seen < 5:
            x, seen = parent_of.get(x, ""), seen + 1
        return x if x in genes else None

    exons = defaultdict(list)
    for par, iv in exons_raw.items():
        g = par if par in genes else to_gene(par)
        if g:
            exons[g].extend(iv)
    exons = {g: merge(iv) for g, iv in exons.items()}
    return genes, tx5, exons


class Blockers:
    """Intervals per chromosome (gene spans, or exons), for overlap queries."""
    def __init__(self, genes, exons=None):
        by = defaultdict(list)
        for g in genes.values():
            if exons is None:
                by[g["chrom"]].append((g["start"], g["end"], g["id"], g["strand"], g["five"]))
            else:
                for s, e in exons.get(g["id"], [(g["start"], g["end"])]):
                    by[g["chrom"]].append((s, e, g["id"], g["strand"], g["five"]))
        self.c = {}
        for chrom, iv in by.items():
            iv.sort()
            starts = [x[0] for x in iv]
            pmax, m = [], 0
            for x in iv:
                m = max(m, x[1])
                pmax.append(m)
            self.c[chrom] = (iv, starts, pmax)

    def overlapping(self, chrom, a, b, exclude):
        iv, starts, pmax = self.c[chrom]
        j = bisect.bisect_right(starts, b) - 1
        while j >= 0 and pmax[j] >= a:
            s, e, gid, strand, five = iv[j]
            if e >= a and gid != exclude:
                yield s, e, gid, strand, five
            j -= 1


def own_introns(g, exons):
    """intron intervals = gene span minus union of its exons."""
    ex = exons.get(g["id"])
    if not ex:
        return []
    out, cur = [], g["start"]
    for s, e in ex:
        if s > cur:
            out.append((cur, s - 1))
        cur = max(cur, e + 1)
    if cur <= g["end"]:
        out.append((cur, g["end"]))
    return out


def upstream_segment(chrom, strand, tss, gid, blk, introns, use_r2):
    """Return (seg_start, seg_end, length, reason, blocker, len_without_r2)."""
    if strand == "+":
        a, b = tss - L_MAX, tss - 1
    else:
        a, b = tss + 1, tss + L_MAX
    a = max(a, 1)
    # R1: other genes. Nearest blocking coordinate on the anchor side.
    cut, who = None, ""
    for s, e, oid, _, _ in blk.overlapping(chrom, a, b, gid):
        edge = e if strand == "+" else s
        edge = min(max(edge, a - 1), b + 1)
        if strand == "+":
            if cut is None or edge > cut:
                cut, who = edge, oid
        else:
            if cut is None or edge < cut:
                cut, who = edge, oid
    r1_a, r1_b = (cut + 1 if (strand == "+" and cut is not None) else a,
                  cut - 1 if (strand == "-" and cut is not None) else b)
    len_r1 = max(0, r1_b - r1_a + 1)
    reason = "neighbour_gene" if cut is not None else ""
    # R2: own introns
    if use_r2:
        for s, e in introns:
            if e < r1_a or s > r1_b:
                continue
            if strand == "+":
                if e + 1 > r1_a:
                    r1_a, reason, who = e + 1, "own_intron", gid
            else:
                if s - 1 < r1_b:
                    r1_b, reason, who = s - 1, "own_intron", gid
    length = max(0, r1_b - r1_a + 1)
    return r1_a, r1_b, length, reason if length < L_MAX else "", who if length < L_MAX else "", len_r1


def inside_other_gene(chrom, tss, gid, blk):
    for s, e, oid, _, _ in blk.overlapping(chrom, tss, tss, gid):
        return oid
    return ""


def head_to_head(chrom, strand, tss, gid, blk, genes):
    a, b = (tss - H2H_WINDOW, tss - 1) if strand == "+" else (tss + 1, tss + H2H_WINDOW)
    for s, e, oid, ostrand, five in blk.overlapping(chrom, max(a, 1), b, gid):
        if ostrand != strand and a <= five <= b:
            return oid
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--cage", action="store_true", help="use FANTOM5 CAGE start sites for genes without EPDnew")
    a = ap.parse_args()
    root = Path(a.root)
    raw, qc, out = root / "datasets" / "raw", root / "datasets" / "qc", root / "datasets" / "cohorts"
    out.mkdir(parents=True, exist_ok=True)

    print("Reading GFF3 ...")
    genes, tx5, exons = load_gff(raw / GFF)
    blk = Blockers(genes)                 # other genes' full spans (rule R1)
    blk_exon = Blockers(genes, exons)     # other genes' exons only (C4: host introns allowed)
    print(f"  {len(genes):,} genes/pseudogenes; exon models for {len(exons):,}")

    # EPDnew promoters per gene (matching done by compare_tss_anchors.py)
    epd = defaultdict(list)
    with open(qc / EPD_MATCH, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["gff_gene_id"] in genes:
                epd[r["gff_gene_id"]].append((r["epd_id"], int(r["epd_tss"]), r["is_main_promoter"] == "True"))

    cage = {}
    if a.cage:
        with open(qc / "cage_anchor_match.csv", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r["gene_id"] in genes and not r["epd_main_tss"]:
                    cage[r["gene_id"]] = (r["cage_peak"], int(r["cage_tss"]))
        print(f"  CAGE start sites for {len(cage):,} genes without EPDnew")

    rows = []
    for gid, g in genes.items():
        chrom, strand = g["chrom"], g["strand"]
        introns = own_introns(g, exons)
        anchors = []
        if gid in epd:
            plist = sorted(epd[gid], key=lambda x: (not x[2], x[0]))
            for i, (eid, tss, main) in enumerate(plist):
                anchors.append((eid, tss, "primary" if i == 0 else "alternative", "EPDnew"))
        elif gid in cage:
            pid, ctss = cage[gid]
            anchors.append((gid + "|cage", ctss, "primary", "CAGE"))
            for t in sorted(set(tx5.get(gid, ())) | {g["five"]}):
                if t != ctss:
                    anchors.append((f"{gid}|tx{t}" if t != g["five"] else gid + "|edge", t, "alternative", "annotation"))
        else:
            anchors.append((gid + "|edge", g["five"], "primary", "annotation"))
            for t in sorted(tx5.get(gid, ())):
                if t != g["five"]:
                    anchors.append((f"{gid}|tx{t}", t, "alternative", "annotation"))
        for aid, tss, atype, src in anchors:
            nested = inside_other_gene(chrom, tss, gid, blk)
            s, e, length, reason, who, len_no_r2 = upstream_segment(
                chrom, strand, tss, gid, blk_exon if nested else blk, introns,
                use_r2=(atype == "alternative"))
            h2h = head_to_head(chrom, strand, tss, gid, blk, genes) if atype == "primary" else ""
            if nested:
                cohort = "C4"
            elif atype == "primary":
                cohort = "C0" if length == 0 else "C1" if length == L_MAX else ("C2b" if h2h else "C2")
            else:
                cohort = "C3" if length >= MIN_ALT else "C3x"
            rows.append(dict(
                anchor_id=aid, gene_id=gid, symbol=g["symbol"], gene_class=g["gclass"],
                biotype=g["biotype"], chrom=chrom, strand=strand, anchor_tss=tss,
                anchor_type=atype, anchor_source=src, cohort=cohort,
                seg_start=s if length else "", seg_end=e if length else "", seg_len=length,
                len_before_R2=len_no_r2, trunc_reason=reason, blocker_id=who,
                nested_in=nested, head_to_head_partner=h2h,
                offset_from_gene_edge=(tss - g["five"]) if strand == "+" else (g["five"] - tss)))

    cols = list(rows[0])
    with open(out / "hpsad_anchors_v7.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    def tab(sel):
        t = defaultdict(Counter)
        for r in sel:
            t[r["cohort"]][r["gene_class"]] += 1
        return {k: dict(v) for k, v in sorted(t.items())}

    prim = [r for r in rows if r["anchor_type"] == "primary"]
    summary = {
        "genes": len(genes),
        "anchors": len(rows),
        "cohort_by_class_all": tab(rows),
        "cohort_by_class_EPDnew_anchored": tab([r for r in rows if r["anchor_source"] == "EPDnew"]),
        "cohort_by_class_annotation_anchored": tab([r for r in rows if r["anchor_source"] == "annotation"]),
        "cohort_by_class_CAGE_anchored": tab([r for r in rows if r["anchor_source"] == "CAGE"]),
        "primary_trunc_reason": dict(Counter(r["trunc_reason"] or "none" for r in prim)),
        "C4_primary_full_length": sum(1 for r in prim if r["cohort"] == "C4" and r["seg_len"] == L_MAX),
        "parameters": {"L_MAX": L_MAX, "H2H_WINDOW": H2H_WINDOW, "MIN_ALT": MIN_ALT},
    }
    with open(out / "hpsad_anchors_v7_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
