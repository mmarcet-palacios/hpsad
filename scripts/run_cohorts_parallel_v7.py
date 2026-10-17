#!/usr/bin/env python3
"""
run_cohorts_parallel_v7.py  -  motif discovery for many cohorts at once

Collects every discovery run (real data + 5 nulls, core and distal zones) for all
cohorts into ONE process pool, largest cohorts first, so all cores stay busy.
Null results cached for the same label AND tag are reused (cache names carry the tag). After all runs
finish, each cohort gets the same outputs as hpsad_motifs_v7.py (FDR per zone and
window class) and is then annotated with hpsad_motif_annotation_v7.py.

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\run_cohorts_parallel_v7.py --jobs 20 --tag _v2
    python scripts\\run_cohorts_parallel_v7.py --jobs 20 --tag _v2 --only C1_transcribed,C4_transcribed
    python scripts\\hpsad_coverage_v7.py --label C1_transcribed   (writes the unexplained id lists)
    python scripts\\run_cohorts_parallel_v7.py --jobs 20 --tag _v2 --only C1_unexpl_cpg,C1_unexpl_noncpg
  stability replicates (independent nulls and composition reference):
    python scripts\\run_cohorts_parallel_v7.py --jobs 20 --nulls 10 --only C1_evidence --tag _v3s1 --seed-base 11000
    (then _v3s2 --seed-base 12000, _v3s3 --seed-base 13000), then scripts\\hpsad_stability_v7.py
  or all replicates in ONE pool (faster: no idle workers between replicates):
    python scripts\\run_cohorts_parallel_v7.py --gpu --jobs 22 --nulls 10 --only C1_evidence --replicates _v5s1:11000,_v5s2:12000,_v5s3:13000
Memory: about 4 GB per running job for C1-sized cohorts (less for small ones).
"""

import argparse
import concurrent.futures as cf
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

if "--gpu" in sys.argv:                     # must be set before the workers import hpsad_motifs_v7
    os.environ["HPSAD_GPU"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402

COHORTS = [  # label, cohort, gene_class, primary_only
    ("C1", "C1", None, False),
    ("C1_protein_coding", "C1", "protein_coding", False),
    ("C1_non_coding", "C1", "non_coding", False),
    ("C1_pseudogene", "C1", "pseudogene", False),
    ("C2b", "C2b", None, False),
    ("C2", "C2", None, False),
    ("C3", "C3", None, False),
    ("C4_primary", "C4", None, True),
    ("C1_transcribed", "C1", "!pseudogene", False),     # C1 without pseudogenes (no -17 scar)
    ("C4_transcribed", "C4", "!pseudogene", True),
    # promoters with no specific element after coverage scan (hpsad_coverage_v7.py), by CpG class
    ("C1_unexpl_cpg", "C1", "!pseudogene", False, "datasets/coverage/ids_C1_transcribed_unexplained_cpg.txt"),
    ("C1_unexpl_noncpg", "C1", "!pseudogene", False, "datasets/coverage/ids_C1_transcribed_unexplained_noncpg.txt"),
    # promoters with initiation evidence only (anchor from EPDnew or FANTOM5 CAGE; needs build_cohorts_v7.py --cage)
    ("C1_evidence", "C1", "!pseudogene", False, {"sources": ["EPDnew", "CAGE"]}),
    ("C4_evidence", "C4", "!pseudogene", True, {"sources": ["EPDnew", "CAGE"]}),
    # -300..+100 around the start site (extract_downstream_v7.py): composite sequences, positions shifted by 100;
    # only the core zone, which now spans depths 1..400 = positions +100..-300
    ("C1ev_down", "C1", "!pseudogene", False,
     {"sources": ["EPDnew", "CAGE"], "segments": "hpsad_segments_v7_down100.tsv.gz", "shift": 100, "zones": {"core": [1, 400]}}),
    ("C1ev_down_tata_pos", "C1", "!pseudogene", False,
     {"sources": ["EPDnew", "CAGE"], "ids_file": "datasets/profiles/ids_C1_evidence_tata_pos.txt",
      "segments": "hpsad_segments_v7_down100.tsv.gz", "shift": 100, "zones": {"core": [1, 400]}}),
    ("C1ev_down_tata_neg_cpg", "C1", "!pseudogene", False,
     {"sources": ["EPDnew", "CAGE"], "ids_file": "datasets/profiles/ids_C1_evidence_tata_neg_cpg.txt",
      "segments": "hpsad_segments_v7_down100.tsv.gz", "shift": 100, "zones": {"core": [1, 400]}}),
    ("C1ev_down_tata_neg_noncpg", "C1", "!pseudogene", False,
     {"sources": ["EPDnew", "CAGE"], "ids_file": "datasets/profiles/ids_C1_evidence_tata_neg_noncpg.txt",
      "segments": "hpsad_segments_v7_down100.tsv.gz", "shift": 100, "zones": {"core": [1, 400]}}),
    # C1_evidence split by TATA status and CpG class (hpsad_syntax_v7.py --write-ids)
    ("C1ev_tata_pos", "C1", "!pseudogene", False,
     {"sources": ["EPDnew", "CAGE"], "ids_file": "datasets/profiles/ids_C1_evidence_tata_pos.txt"}),
    ("C1ev_tata_neg_cpg", "C1", "!pseudogene", False,
     {"sources": ["EPDnew", "CAGE"], "ids_file": "datasets/profiles/ids_C1_evidence_tata_neg_cpg.txt"}),
    ("C1ev_tata_neg_noncpg", "C1", "!pseudogene", False,
     {"sources": ["EPDnew", "CAGE"], "ids_file": "datasets/profiles/ids_C1_evidence_tata_neg_noncpg.txt"}),
]
# Figure 7 cohort comparison with the current code, on -300..+100 sequences (labels end in _d):
# core zone = -300..+100, distal zone = -301..-1900. Family weights from hpsad_families_core_v7.py (copies shared
# only within -300..+100, e.g. multicopy small-RNA gene bodies, count once). Gene classes: small-RNA genes
# (tRNA, snRNA, snoRNA, miRNA ...) are class 'small_rna' (hpsad_motifs_v7.SMALL_RNA), no longer part of 'non_coding'.
_DOWN = {"segments": "hpsad_segments_v7_down100.tsv.gz", "families": "hpsad_families_v7_down100.tsv", "shift": 100,
         "zones": {"core": [1, 400], "distal": [401, 2000]}}
COHORTS += [(lab + "_d", coh, gc, prim, dict(_DOWN)) for lab, coh, gc, prim in [
    ("C1", "C1", None, False), ("C1_protein_coding", "C1", "protein_coding", False), ("C1_non_coding", "C1", "non_coding", False),
    ("C1_pseudogene", "C1", "pseudogene", False), ("C2", "C2", None, False), ("C2b", "C2b", None, False), ("C3", "C3", None, False),
    ("C4_primary", "C4", None, True), ("C1_transcribed", "C1", "!pseudogene", False), ("C4_transcribed", "C4", "!pseudogene", True),
    ("C1_small_rna", "C1", "small_rna", False), ("C4_small_rna", "C4", "small_rna", True),
    ("C1_pol2", "C1", "!pseudogene,small_rna", False), ("C4_pol2", "C4", "!pseudogene,small_rna", True)]]


def position(d, shift):
    """start depth -> position relative to the start site (no 0): upstream -1, -2, ...; downstream +1, +2, ..."""
    return -(d - shift) if d > shift else shift + 1 - d


def rows_of(motifs, keep, zone, shift=0):
    rows = []
    for m in motifs:
        known = [nm for nm, cons, (lo, hi) in M.KNOWN
                 if lo - 5 <= m["d_start"] - shift <= hi + 15 and M.core_match(m["consensus"], cons)]
        p = position(m["d_start"], shift)
        rows.append({"zone": zone, "rank": m["rank"], "pass_fdr": m["rank"] in keep, "consensus": m["consensus"],
                     "start_depth": m["d_start"], "position": f"{p:+d}", "length": m["length"], "window": m["width"],
                     "O": round(m["O"], 1), "E": round(m["E"], 1), "ratio": round(m["O"] / max(m["E"], 1e-9), 2),
                     "score": round(m["score"], 1), "instances": m["n_inst"],
                     "weighted_instances": round(m["w_inst"], 1), "seed": m["seed"][1], "known": ",".join(known),
                     "pwm_ACGT": ";".join(",".join(f"{x:.3f}" for x in row) for row in m["prob"])})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--jobs", type=int, default=20)
    ap.add_argument("--workers", type=int, default=16, help="processes for shuffling")
    ap.add_argument("--nulls", type=int, default=5)
    ap.add_argument("--max-motifs", type=int, default=M.MAX_MOTIFS)
    ap.add_argument("--tag", default="_v2")
    ap.add_argument("--only", default="", help="comma-separated labels to run (default: all)")
    ap.add_argument("--gpu", action="store_true", help="word scan on the GPU (CuPy); use fewer --jobs, e.g. 8")
    ap.add_argument("--seed-base", type=int, default=2000,
                    help="null shuffles use seeds seed-base+r, the composition reference seed-base+999; "
                         "use a different value (and tag) for each stability replicate")
    ap.add_argument("--replicates", default="",
                    help="several replicates in one pool: tag:seed-base,tag:seed-base,... (overrides --tag and --seed-base)")
    a = ap.parse_args()
    reps = [(t, int(b)) for t, b in (x.split(":") for x in a.replicates.split(","))] if a.replicates else [(a.tag, a.seed_base)]
    root = Path(a.root)
    out = root / "datasets" / "motifs"
    (out / "nulls").mkdir(parents=True, exist_ok=True)
    T0 = time.time()
    todo = [c for c in COHORTS if not a.only or c[0] in a.only.split(",")]

    ctx, jobs = {}, []
    for (label0, cohort, gclass, prim, *extra), (tag, seed_base) in [(c, r) for c in todo for r in reps]:
        opt = extra[0] if extra else {}
        opt = {"ids_file": opt} if isinstance(opt, str) else opt
        ids, seqs, w, cls = M.load_cohort(root, cohort, gclass, prim, ids_file=opt.get("ids_file"),
                                          sources=opt.get("sources"),
                                          segments=opt.get("segments", "hpsad_segments_v7.tsv.gz"),
                                          families=opt.get("families", "hpsad_families_v7_down100.tsv"))
        zones = {z: tuple(v) for z, v in opt.get("zones", M.ZONES).items()}      # depth ranges (shifted cohorts differ)
        label = (label0, tag)
        if len(seqs) < 200:
            print(f"{label}: only {len(seqs)} promoters - skipped")
            continue
        ref = M.null_cohort(seqs, seed_base + 999, out / "nulls" / f"null50_{label0}{tag}_reference.txt.gz", a.workers)
        Sref = M.encode(ref)
        nulls = {z: [None] * a.nulls for z in zones}
        for z, zr in zones.items():
            for r in range(a.nulls):
                cache = out / "nulls" / f"motifs_null_{label0}{tag}_{z}_rep{r}_m{a.max_motifs}.json"
                if cache.exists():
                    nulls[z][r] = json.load(open(cache))
                    continue
                ns = M.null_cohort(seqs, seed_base + r, out / "nulls" / f"null50_{label0}{tag}_rep{r}.txt.gz", a.workers)
                jobs.append((len(seqs), (f"{label0}{tag}|{z}:null{r}", M.encode(ns), w, Sref, a.max_motifs, zr)))
            jobs.append((len(seqs), (f"{label0}{tag}|{z}:real", M.encode(seqs), w, Sref, a.max_motifs, zr)))
        ctx[label] = {"n": len(seqs), "w": w, "nulls": nulls, "zones": zones, "shift": opt.get("shift", 0)}
        print(f"{label0}{tag}: {len(seqs):,} promoters, {w.sum():,.0f} independent "
              f"({sum(1 for _, j in jobs if j[0].startswith(label0 + tag + '|'))} runs queued)")
    jobs.sort(key=lambda x: -x[0])                                  # largest first
    print(f"\n{len(jobs)} discovery runs, {a.jobs} at a time ({time.time() - T0:.0f} s to prepare)"
          f"{'; word scan on GPU' if a.gpu else ''}")

    results = {}
    with cf.ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(M._discover_job, j): j[0] for _, j in jobs}
        for k, f in enumerate(cf.as_completed(futs), 1):
            tag, ms = f.result()
            results[tag] = ms
            print(f"  [{k}/{len(jobs)}] finished {tag}: {len(ms)} motifs ({(time.time() - T0) / 60:.0f} min)")

    for (label, tag), c in ctx.items():
        for z in c["zones"]:
            for r in range(a.nulls):
                if c["nulls"][z][r] is None:
                    slim = [{k: v for k, v in m.items() if k in ("score", "width", "d_start", "consensus", "O", "E", "rank")}
                            for m in results[f"{label}{tag}|{z}:null{r}"]]
                    json.dump(slim, open(out / "nulls" / f"motifs_null_{label}{tag}_{z}_rep{r}_m{a.max_motifs}.json", "w"))
                    c["nulls"][z][r] = slim
        real_rows = []
        for z in c["zones"]:
            real = results[f"{label}{tag}|{z}:real"]
            real_rows += rows_of(real, M.fdr_select(real, c["nulls"][z]), z, c["shift"])
        path = out / f"hpsad_motifs_{label}_zones{tag}.tsv"
        with open(path, "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(real_rows[0]) if real_rows else ["zone"], delimiter="\t")
            wr.writeheader()
            wr.writerows(real_rows)
        summary = {"cohort": label, "promoters": c["n"], "independent": round(float(c["w"].sum()), 1),
                   "nulls": a.nulls, "shift": c["shift"], "zones": c["zones"],
                   "null_max_scores": {z: [max([m["score"] for m in nm] or [0]) for nm in c["nulls"][z]] for z in c["zones"]},
                   "motifs": real_rows}
        json.dump(summary, open(out / f"hpsad_motifs_{label}_zones{tag}_summary.json", "w"), indent=1, default=float)
        npass = {z: sum(r["pass_fdr"] for r in real_rows if r["zone"] == z) for z in c["zones"]}
        print(f"\n{label}{tag}: pass FDR - " + ", ".join(f"{z} {n}" for z, n in npass.items()))
        for r in real_rows:
            if r["pass_fdr"]:
                print(f"  {r['zone']:<6} M{r['rank']:02d} {r['consensus']:<18} {r['position']:<6} w={r['window']:<3} "
                      f"x{r['ratio']:<6} score {r['score']:<7} n={r['instances']:<6} {r['known']}")
        subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "hpsad_motif_annotation_v7.py"),
                        "--root", str(root), "--cohort", label, "--tag", tag], check=False)
    print(f"\nAll done in {(time.time() - T0) / 60:.0f} min")


if __name__ == "__main__":
    main()
