#!/usr/bin/env python3
"""
hpsad_tsd_v7.py  -  target-site duplications (TSDs) flanking genes: direct test of the -17 scar

A gene copy pasted into the genome by LINE-1 (a processed pseudogene) is flanked by a short
direct repeat: the target site duplication. One copy sits just upstream of the 5' edge, the
other just downstream of the 3' end (after the poly(A) tail). This script looks for that repeat.

For each gene whose primary anchor is in C1, C2 or C4 (gene edges from the GFF, gene orientation):
    up   = 40 nt upstream of the 5' edge (-40..-1)
    down = 120 nt downstream of the 3' end (+1..+120), long enough to pass a poly(A) tail
TSD = longest exact match (>= MIN_TSD nt) between up[-30:] and down, not low-complexity
      (no base > 75% of the match). Null: the same search with up and down from DIFFERENT
      genes of the same group (5 random pairings) - the chance rate for that composition.
Also scored: the L1 target-site word TTAAAR/TTRAAA starting at -25..-10 ("scar").

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_tsd_v7.py
Outputs: datasets\\qc\\tsd_v7_summary.tsv  (per cohort x gene class)
         datasets\\qc\\tsd_v7_genes.tsv    (per gene: TSD length, position, sequence, scar)
"""
import argparse
import collections
import csv
import gzip
import random
import re
from pathlib import Path

FASTA = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
GFF = "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
UP, UP_SEARCH, DOWN = 40, 30, 120
MIN_TSD = 12
NULL_REPS = 5
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SCAR = re.compile(r"(?=(?:TT[AG]AAA|TTAAA[AG]))")


def rc(s):
    return s.translate(COMP)[::-1]


def stream_fasta(path, wanted):
    name, parts = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            elif name in wanted:
                parts.append(line.rstrip())
        if name in wanted:
            yield name, "".join(parts)


def gene_coords(path, wanted):
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            m = re.search(r"ID=([^;\n]+)", f[8])
            if m and m.group(1) in wanted:
                out[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
    return out


def tsd(up, down, min_len=MIN_TSD):
    """longest exact, non-low-complexity match between up[-UP_SEARCH:] and down.
    Returns (length, depth of 5' copy start (positive, from the 5' edge), sequence, offset in down)."""
    a = up[-UP_SEARCH:]
    for k in range(len(a), min_len - 1, -1):
        for i in range(len(a) - k + 1):
            w = a[i:i + k]
            if "N" in w or max(w.count(b) for b in "ACGT") > 0.75 * k:
                continue
            j = down.find(w)
            if j >= 0:
                return k, len(a) - i, w, j
    return 0, None, "", None


def has_scar(up):
    L = len(up)
    return any(10 <= L - m.start() <= 25 for m in SCAR.finditer(up))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--gff", default=None)
    ap.add_argument("--anchors", default=None)
    a = ap.parse_args()
    root = Path(a.root)
    raw, qc = root / "datasets" / "raw", root / "datasets" / "qc"
    qc.mkdir(parents=True, exist_ok=True)

    genes = {}
    with open(a.anchors or root / "datasets" / "cohorts" / "hpsad_anchors_v7.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["anchor_type"] == "primary" and r["cohort"] in ("C1", "C2", "C4"):
                genes[r["gene_id"]] = dict(cohort=r["cohort"], cls=r["gene_class"], symbol=r["symbol"])
    print(f"{len(genes):,} genes with a primary anchor in C1, C2 or C4; reading GFF ...")
    coords = gene_coords(a.gff or raw / GFF, set(genes))
    by_chrom = collections.defaultdict(list)
    for g, c in coords.items():
        by_chrom[c[0]].append(g)
    print(f"{len(coords):,} located; reading genome ...")
    for chrom, seq in stream_fasta(a.fasta or raw / FASTA, set(by_chrom)):
        for g in by_chrom[chrom]:
            _, s, e, strand = coords[g]
            if strand == "+":
                up, down = seq[max(0, s - 1 - UP):s - 1], seq[e:e + DOWN]
            else:
                up, down = rc(seq[e:e + UP]), rc(seq[max(0, s - 1 - DOWN):s - 1])
            genes[g]["up"], genes[g]["down"] = up.upper(), down.upper()
        print(f"  {chrom} done")

    groups = collections.defaultdict(list)
    for g, d in genes.items():
        if "up" not in d or len(d["up"]) < UP or len(d["down"]) < DOWN:
            continue
        d["tsd"] = tsd(d["up"], d["down"])
        d["scar"] = has_scar(d["up"])
        groups[(d["cohort"], d["cls"])].append(g)
        if d["cls"] == "pseudogene" and re.match(r"RP[LS]\d", d["symbol"]):
            groups[(d["cohort"], "pseudogene_RPL_RPS")].append(g)

    rng = random.Random(7)
    rows = []
    for (c, cl), ids in sorted(groups.items()):
        n = len(ids)
        hit = [genes[g]["tsd"][0] >= MIN_TSD for g in ids]
        null = 0.0
        for _ in range(NULL_REPS):
            perm = ids[:]
            rng.shuffle(perm)
            null += sum(tsd(genes[x]["up"], genes[y]["down"])[0] >= MIN_TSD
                        for x, y in zip(ids, perm[1:] + perm[:1])) / n
        null /= NULL_REPS
        scar = [genes[g]["scar"] for g in ids]
        tsd_scar = sum(h and s for h, s in zip(hit, scar)) / max(sum(scar), 1)
        tsd_noscar = sum(h and not s for h, s in zip(hit, scar)) / max(n - sum(scar), 1)
        depths = sorted(genes[g]["tsd"][1] for g, h in zip(ids, hit) if h)
        lens = sorted(genes[g]["tsd"][0] for g, h in zip(ids, hit) if h)
        rows.append([c, cl, n, sum(hit) / n, null, (sum(hit) / n) / null if null > 0 else "NA", sum(scar) / n,
                     tsd_scar, tsd_noscar,
                     lens[len(lens) // 2] if lens else "", depths[len(depths) // 2] if depths else ""])

    hdr = ["cohort", "group", "n", "tsd_frac", "tsd_null", "tsd_fold", "scar_frac",
           "tsd_frac_if_scar", "tsd_frac_if_no_scar", "median_tsd_len", "median_tsd_start_depth"]
    with open(qc / "tsd_v7_summary.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        for r in rows:
            w.writerow([f"{x:.4f}" if isinstance(x, float) else x for x in r])
    with open(qc / "tsd_v7_genes.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "symbol", "cohort", "gene_class", "tsd_len", "tsd_start_depth",
                     "tsd_seq", "tsd_offset_downstream", "scar_-25_-10"])
        for g, d in genes.items():
            if "tsd" in d:
                k, dep, seq, off = d["tsd"]
                w.writerow([g, d["symbol"], d["cohort"], d["cls"], k, dep or "", seq,
                            "" if off is None else off, int(d["scar"])])

    print(f"\nTSD = exact repeat >= {MIN_TSD} nt shared by -30..-1 and +1..+{DOWN}; null = mismatched gene pairs")
    print(f"{'cohort':<7}{'group':<20}{'n':>7}{'TSD':>7}{'null':>7}{'fold':>6}{'scar':>7}"
          f"{'TSD|scar':>10}{'TSD|no':>8}{'len':>5}{'depth':>6}")
    for r in rows:
        print(f"{r[0]:<7}{r[1]:<20}{r[2]:>7}{r[3]:>7.3f}{r[4]:>7.3f}{(f'{r[5]:.1f}' if r[5] != 'NA' else 'NA'):>6}{r[6]:>7.3f}"
              f"{r[7]:>10.3f}{r[8]:>8.3f}{str(r[9]):>5}{str(r[10]):>6}")
    print(f"\nwritten {qc / 'tsd_v7_summary.tsv'} and tsd_v7_genes.tsv")


if __name__ == "__main__":
    main()
