#!/usr/bin/env python3
"""
hpsad_l1_footprint_v7.py  -  is the AT-rich motif at -17 a retrotransposition footprint?

For every promoter segment (>= 400 nt, no N) of C1, C2, C2b, C3 and C4 (primary anchors),
records whether the L1 target-site word TTAAAR / TTRAAA starts at -25..-10 (16 positions)
and, as background, in 16-nt windows starting at -200, -250, -300, -350 (averaged).
TATAWA is scored the same way at -35..-25 (11-nt windows).
Groups: cohort x gene class, plus RPL/RPS pseudogenes (classic processed pseudogenes).

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_l1_footprint_v7.py
Output: datasets\\qc\\l1_footprint_v7.tsv
"""
import argparse
import collections
import csv
import gzip
import re
from pathlib import Path

L1 = re.compile(r"(?=(?:TT[AG]AAA|TTAAA[AG]))")
TATA = re.compile(r"(?=TATA[AT]A)")
BG = (200, 250, 300, 350)


def depths(s, pat):
    L = len(s)
    return {L - m.start() for m in pat.finditer(s)}


def any_in(ds, lo, hi):
    return any(lo <= d <= hi for d in ds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    root = Path(a.root)
    csv.field_size_limit(10 ** 8)
    acc = collections.defaultdict(lambda: [0, 0, 0.0, 0, 0.0, 0])  # n, l1, l1_bg, tata, tata_bg, both
    with gzip.open(root / "datasets" / "cohorts" / "hpsad_segments_v7.tsv.gz", "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = r["cohort"]
            if c not in ("C1", "C2", "C2b", "C3", "C4") or r["has_N"] != "0":
                continue
            if c == "C4" and r["anchor_type"] != "primary":
                continue
            s = r["sequence_5to3"].upper()
            if len(s) < 400:
                continue
            dl, dt = depths(s, L1), depths(s, TATA)
            l1 = any_in(dl, 10, 25)
            ta = any_in(dt, 25, 35)
            keys = [(c, r["gene_class"])]
            if r["gene_class"] == "pseudogene" and re.match(r"RP[LS]\d", r["symbol"]):
                keys.append((c, "pseudogene_RPL_RPS"))
            for k in keys:
                v = acc[k]
                v[0] += 1
                v[1] += l1
                v[2] += sum(any_in(dl, o, o + 15) for o in BG) / len(BG)
                v[3] += ta
                v[4] += sum(any_in(dt, o, o + 10) for o in BG) / len(BG)
                v[5] += l1 and ta
    out = root / "datasets" / "qc" / "l1_footprint_v7.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cohort", "group", "n", "l1_frac_-25_-10", "l1_frac_bg", "l1_fold",
                    "tata_frac_-35_-25", "tata_frac_bg", "both_frac"])
        print(f"{'cohort':<6}{'group':<22}{'n':>7}{'L1 -25..-10':>12}{'bg':>7}{'fold':>6}{'TATA':>8}{'bg':>7}{'both':>7}")
        for (c, g), (n, l1, l1b, ta, tab, both) in sorted(acc.items()):
            row = [c, g, n, l1 / n, l1b / n, (l1 / n) / max(l1b / n, 1e-9), ta / n, tab / n, both / n]
            w.writerow([f"{x:.4f}" if isinstance(x, float) else x for x in row])
            print(f"{c:<6}{g:<22}{n:>7}{row[3]:>12.3f}{row[4]:>7.3f}{row[5]:>6.1f}{row[6]:>8.3f}{row[7]:>7.3f}{row[8]:>7.3f}")
    print(f"\nwritten {out}")


if __name__ == "__main__":
    main()
