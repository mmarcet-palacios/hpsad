#!/usr/bin/env python3
"""
extract_downstream_v7.py  -  add the first D nt downstream of each start site (+1..+D)

Reads datasets/cohorts/hpsad_segments_v7.tsv.gz (upstream segments, last base = -1) and the
GRCh38.p14 FASTA, and writes a composite cohort file in which every sequence is
    upstream segment  +  downstream +1..+D   (transcriptional orientation, FASTA case kept)
so the discovery code (which reads the last L = 2000 nt) sees -(L-D)..-1 followed by +1..+D.
Position convention in the composite: +1 is the start-site base itself (anchor_tss); there is
no position 0. A start depth d in the composite (start index j = L - d) corresponds to
    position = -(d - D)   for d > D   (upstream)
    position = D + 1 - d  for d <= D  (downstream)

  datasets/cohorts/hpsad_segments_v7_down<D>.tsv.gz   same columns as the segments file, with
      sequence_5to3 / sequence_softmasked replaced by the composite, has_N updated, plus
      down_seq (the +1..+D part) and shift = D
  datasets/cohorts/hpsad_segments_v7_down<D>_summary.json

Self-check: for every segment the 10 nt (or the whole segment, if shorter) just upstream of anchor_tss are re-read from the genome
and must equal the last 10 nt of the stored upstream segment (catches strand / off-by-one errors).

Usage (from C:\\Promoter\\HPSAD_v7):   python scripts\\extract_downstream_v7.py            (D = 100)
                                         python scripts\\extract_downstream_v7.py --down 50
"""
import argparse
import csv
import gzip
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

FASTA = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
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


def downstream(seq, tss, strand, D):
    """+1..+D in transcriptional orientation (tss is 1-based); also the 10 nt at -10..-1."""
    if strand == "+":
        down = seq[tss - 1:tss - 1 + D]
        up10 = seq[tss - 11:tss - 1]
    else:
        down = revcomp(seq[max(0, tss - D):tss])
        up10 = revcomp(seq[tss:tss + 10])
    return down, up10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--down", type=int, default=100)
    a = ap.parse_args()
    root = Path(a.root)
    raw, coh = root / "datasets" / "raw", root / "datasets" / "cohorts"
    D = a.down
    csv.field_size_limit(10 ** 8)
    by_chrom = defaultdict(list)
    with gzip.open(coh / "hpsad_segments_v7.tsv.gz", "rt", encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        for r in rd:
            by_chrom[r["chrom"]].append(r)
    total = sum(len(v) for v in by_chrom.values())
    print(f"{total:,} segments on {len(by_chrom)} chromosomes; adding +1..+{D}")
    out_cols = cols + ["down_seq", "shift"]
    st, mism, short, t0, done = Counter(), [], 0, time.time(), 0
    out = coh / f"hpsad_segments_v7_down{D}.tsv.gz"
    with gzip.open(out, "wt", encoding="utf-8", newline="") as gz:
        w = csv.DictWriter(gz, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        for chrom, seq in stream_fasta(raw / FASTA, set(by_chrom)):
            for r in by_chrom[chrom]:
                tss = int(r["anchor_tss"])
                down, up10 = downstream(seq, tss, r["strand"], D)
                k = min(10, len(r["sequence_5to3"]))       # segments shorter than 10 nt: compare what exists
                if up10[10 - k:].upper() != r["sequence_5to3"][-k:].upper():
                    mism.append(r["anchor_id"])
                if len(down) < D:                         # chromosome end: pad with N (excluded below)
                    short += 1
                    down = down + "N" * (D - len(down))
                has_n = int(r["has_N"] == "1" or any(ch not in "ACGT" for ch in down.upper()))
                r.update(sequence_5to3=r["sequence_5to3"] + down.upper(),
                         sequence_softmasked=r["sequence_softmasked"] + down,
                         has_N=has_n, down_seq=down, shift=D)
                w.writerow(r)
                st[r["cohort"]] += 1
                st["downstream_with_N"] += any(ch not in "ACGT" for ch in down.upper())
            done += len(by_chrom[chrom])
            print(f"  {chrom}: {len(by_chrom[chrom]):,}  ({done:,}/{total:,}, {time.time() - t0:.0f} s)")
            del seq
    summary = {"segments": done, "shift": D, "by_cohort": dict(st), "chromosome_end_padded": short,
               "upstream_check_mismatches": len(mism), "first_mismatches": mism[:20]}
    json.dump(summary, open(coh / f"hpsad_segments_v7_down{D}_summary.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))
    if mism:
        raise SystemExit(f"{len(mism)} segments failed the upstream check - do not use {out.name}")


if __name__ == "__main__":
    main()
