#!/usr/bin/env python3
"""
extract_cohort_sequences_v7.py  -  HPSAD fresh start, step 2b: sequences

Reads datasets/cohorts/hpsad_anchors_v7.tsv (from build_cohorts_v7.py) and the
GRCh38.p14 FASTA, and writes one sequence per anchor with a segment:

  datasets/cohorts/hpsad_segments_v7.tsv.gz
      all anchor columns, plus:
      sequence_5to3   upper case, transcriptional orientation; the LAST base is
                      position -1 (immediately upstream of the start site)
      sequence_softmasked  same sequence with the FASTA's case kept
                      (lower case = repeat-masked by NCBI)
      has_N           1 if any non-ACGT base (excluded from discovery)
      repeat_frac     fraction of bases soft-masked (lower case) in the FASTA
      gc_frac, cpg_oe GC fraction and CpG observed/expected
  datasets/cohorts/hpsad_segments_v7_summary.json
      counts per cohort x gene_class, N exclusions, positive-control check

Streams the FASTA one chromosome at a time (peak memory ~1 GB).
Usage (from the project root):  python scripts\\extract_cohort_sequences_v7.py
"""

import argparse
import csv
import gzip
import json
import time
from collections import defaultdict, Counter
from pathlib import Path

FASTA = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CONTROLS = {"HBB_1", "HBD_1", "HBG1_1", "HBG2_1"}   # TATA-like box expected near -30


def revcomp(s):
    return s.translate(COMP)[::-1]


def stream_fasta(path, wanted):
    """yield (accession, sequence) for wanted accessions only."""
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


def composition(seq_upper):
    n = len(seq_upper)
    g, c = seq_upper.count("G"), seq_upper.count("C")
    cg = seq_upper.count("CG")
    oe = (cg * n / (c * g)) if c and g else 0.0
    return round((g + c) / n, 4), round(oe, 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    root = Path(a.root)
    raw, coh = root / "datasets" / "raw", root / "datasets" / "cohorts"

    by_chrom = defaultdict(list)
    with open(coh / "hpsad_anchors_v7.tsv", encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        for r in rd:
            if int(r["seg_len"]) > 0:
                by_chrom[r["chrom"]].append(r)
    total = sum(len(v) for v in by_chrom.values())
    print(f"{total:,} segments on {len(by_chrom)} chromosomes")

    out_cols = cols + ["sequence_5to3", "sequence_softmasked", "has_N", "repeat_frac", "gc_frac", "cpg_oe"]
    stats = defaultdict(Counter)
    controls, done, t0 = {}, 0, time.time()
    with gzip.open(coh / "hpsad_segments_v7.tsv.gz", "wt", encoding="utf-8", newline="") as gz:
        w = csv.DictWriter(gz, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        for chrom, seq in stream_fasta(raw / FASTA, set(by_chrom)):
            for r in by_chrom[chrom]:
                s, e = int(r["seg_start"]), int(r["seg_end"])
                frag = seq[s - 1:e]
                if r["strand"] == "-":
                    frag = revcomp(frag)
                if len(frag) != int(r["seg_len"]):
                    raise SystemExit(f"Length mismatch for {r['anchor_id']}: {len(frag)} vs {r['seg_len']}")
                up = frag.upper()
                has_n = int(any(ch not in "ACGT" for ch in up))
                gc, oe = composition(up)
                r.update(sequence_5to3=up, sequence_softmasked=frag, has_N=has_n,
                         repeat_frac=round(sum(ch.islower() for ch in frag) / len(frag), 4),
                         gc_frac=gc, cpg_oe=oe)
                w.writerow(r)
                key = f'{r["cohort"]}|{r["gene_class"]}'
                stats[key]["segments"] += 1
                stats[key]["with_N"] += has_n
                stats[key]["full_2000"] += int(r["seg_len"]) == 2000
                if r["anchor_id"] in CONTROLS:
                    controls[r["anchor_id"]] = up[-40:-19]      # positions -40 .. -20
            done += len(by_chrom[chrom])
            print(f"  {chrom}: {len(by_chrom[chrom]):,} segments  ({done:,}/{total:,}, {time.time() - t0:.0f} s)")
            del seq

    summary = {
        "segments_written": done,
        "by_cohort_and_class": {k: dict(v) for k, v in sorted(stats.items())},
        "N_excluded_total": sum(v["with_N"] for v in stats.values()),
        "positive_controls_-40_to_-20": controls,
        "control_note": "HBB/HBD should contain CATAAAA and HBG1/HBG2 AATAAAA near -31",
    }
    with open(coh / "hpsad_segments_v7_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps({k: summary[k] for k in ("segments_written", "N_excluded_total",
                                              "positive_controls_-40_to_-20")}, indent=2))


if __name__ == "__main__":
    main()
