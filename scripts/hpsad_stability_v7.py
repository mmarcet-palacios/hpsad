#!/usr/bin/env python3
"""
hpsad_stability_v7.py  -  which motifs pass FDR reproducibly across independent replicates?

Each replicate is the same cohort run with its own null shuffles and composition reference
(run_cohorts_parallel_v7.py --seed-base ... --tag ...). Motifs from all replicates are matched:
same zone, start positions within --max-shift nt, and PWMs that agree (mean column Pearson
correlation >= --min-pcc over >= 5 overlapping columns, offsets -6..+6, sense strand).
For each matched motif the table gives, per replicate: PASS (passes FDR), found (discovered but
below the FDR cut, with its score) or absent.

  stable      passes FDR in every replicate
  variable    passes in at least one, not all
  (motifs that never pass are not listed)

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_stability_v7.py --label C1_evidence --tags _v3s1,_v3s2,_v3s3
Output: datasets/motifs/hpsad_stability_<label>.tsv
"""
import argparse
import csv
from pathlib import Path

import numpy as np


def load(path):
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        out.append(dict(zone=r["zone"], rank=int(r["rank"]), passed=r["pass_fdr"] == "True", cons=r["consensus"],
                        d=int(r["start_depth"]), oe=float(r["ratio"]), score=float(r["score"]), n=int(r["instances"]),
                        known=r.get("known", ""),
                        pwm=np.array([[float(x) for x in c.split(",")] for c in r["pwm_ACGT"].split(";")])))
    return out


def informative(p):
    return 2 + np.sum(np.where(p > 0, p * np.log2(np.clip(p, 1e-9, 1)), 0), axis=1) >= 0.5


def similar(a, b, max_shift, min_pcc):
    if a["zone"] != b["zone"] or abs(a["d"] - b["d"]) > max_shift + 6:
        return False
    best = -1.0
    for off in range(-6, 7):                      # b column j aligned to a column j + off
        i0, j0 = max(0, off), max(0, -off)
        k = min(len(a["pwm"]) - i0, len(b["pwm"]) - j0)
        if k < 5:
            continue
        A, B = a["pwm"][i0:i0 + k], b["pwm"][j0:j0 + k]
        keep = informative(A) | informative(B)
        if keep.sum() < 5:
            continue
        pcc = [np.corrcoef(x, y)[0, 1] for x, y in zip(A[keep], B[keep]) if x.std() > 0 and y.std() > 0]
        if len(pcc) >= 5:
            best = max(best, float(np.mean(pcc)))
    return best >= min_pcc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--label", default="C1_evidence")
    ap.add_argument("--tags", default="_v3s1,_v3s2,_v3s3")
    ap.add_argument("--max-shift", type=int, default=10)
    ap.add_argument("--min-pcc", type=float, default=0.85)
    a = ap.parse_args()
    mdir = Path(a.root) / "datasets" / "motifs"
    tags = a.tags.split(",")
    reps = {t: load(mdir / f"hpsad_motifs_{a.label}_zones{t}.tsv") for t in tags}

    clusters = []                                  # each: {tag: motif}
    for t in tags:
        for m in sorted(reps[t], key=lambda x: -x["score"]):
            for c in clusters:
                if t not in c and any(similar(m, o, a.max_shift, a.min_pcc) for o in c.values()):
                    c[t] = m
                    break
            else:
                clusters.append({t: m})
    rows = []
    for c in clusters:
        npass = sum(1 for t in tags if t in c and c[t]["passed"])
        if npass == 0:
            continue
        rep = max(c.values(), key=lambda m: m["score"])
        cells = []
        for t in tags:
            m = c.get(t)
            cells.append("absent" if m is None else (f"PASS {m['score']:.0f}" if m["passed"] else f"found {m['score']:.0f}"))
        rows.append(dict(status="stable" if npass == len(tags) else "variable", passes=f"{npass}/{len(tags)}",
                         zone=rep["zone"], consensus=rep["cons"],
                         start_depths=",".join(str(c[t]["d"]) if t in c else "-" for t in tags),
                         oe_mean=round(float(np.mean([m["oe"] for m in c.values()])), 2),
                         instances_mean=int(np.mean([m["n"] for m in c.values()])), known=rep["known"],
                         **{t.strip("_"): x for t, x in zip(tags, cells)}))
    rows.sort(key=lambda r: (r["status"] != "stable", r["zone"], int(r["start_depths"].split(",")[0]) if r["start_depths"].split(",")[0] != "-" else 9999))
    out = mdir / f"hpsad_stability_{a.label}.tsv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{a.label}: replicates {', '.join(tags)}; passing per replicate: "
          + ", ".join(f"{t} {sum(m['passed'] for m in reps[t])}" for t in tags))
    print(f"{'status':<9}{'pass':<6}{'zone':<7}{'consensus':<18}{'starts':<16}{'O/E':>6}{'n':>6}  " + "  ".join(f"{t:<10}" for t in tags))
    for r in rows:
        print(f"{r['status']:<9}{r['passes']:<6}{r['zone']:<7}{r['consensus']:<18}{r['start_depths']:<16}{r['oe_mean']:>6}"
              f"{r['instances_mean']:>6}  " + "  ".join(f"{r[t.strip('_')]:<10}" for t in tags) + f"  {r['known']}")
    print(f"\nstable {sum(r['status'] == 'stable' for r in rows)}, variable {sum(r['status'] == 'variable' for r in rows)}; written {out}")


if __name__ == "__main__":
    main()
