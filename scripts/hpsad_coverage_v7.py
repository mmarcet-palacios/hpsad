#!/usr/bin/env python3
"""
hpsad_coverage_v7.py  -  how many promoters do the discovered elements explain?

For one motif run (e.g. C1_transcribed), every motif passing FDR is scanned back over the
cohort at its own position (start depth +/- window/2 + 3 nt, sense strand, repeat-masked
bases excluded). Its score threshold is calibrated so the number of promoters scoring above
it equals the motif's reported instance count. Motifs are merged into their elements
(hpsad_elements_<label>_v2.tsv). The same scan on a dinucleotide shuffle of each promoter
(50-nt blocks) gives the chance rate.

  specific motif    = a motif with O/E >= --oe-specific (default 2.5); only specific motifs
                      "explain" a promoter (an element counts through its specific motifs only;
                      low-specificity GC-rich motifs are reported but never explain)
  copy artefact     = a motif with O/E > 20 and fewer than 60 carriers: skipped
  explained         = carries >= 1 specific motif
  coverage          = fraction explained (real) minus fraction explained (shuffled)
  CpG class         = -500..-1: GC >= 0.55 and CpG obs/exp >= 0.6 -> "cpg", else "noncpg"

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_coverage_v7.py --label C1_transcribed --cohort C1 --exclude-class pseudogene
    python scripts\\hpsad_coverage_v7.py --label C1_evidence --tag _v3 --sources EPDnew,CAGE
Outputs (datasets\\coverage):
    coverage_<label>.tsv                       one row per promoter
    coverage_<label>_summary.json              per element real vs shuffled, coverage by CpG class
    ids_<label>_unexplained_cpg.txt / _noncpg.txt   anchor ids for the next discovery round
"""
import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402

IDX = {c: i for i, c in enumerate("ACGT")}


def cpg_class(seq):
    u = seq.upper()[-500:]
    n = len(u)
    c, g, cg = u.count("C"), u.count("G"), u.count("CG")
    oe = cg * n / (c * g) if c and g else 0.0
    return "cpg" if (c + g) / n >= 0.55 and oe >= 0.6 else "noncpg"


def encode(seqs, L):
    A = np.full((len(seqs), L), 4, np.int8)
    for k, s in enumerate(seqs):
        s = s[-L:]
        off = L - len(s)
        A[k, off:] = [IDX.get(ch, 4) for ch in s]         # lower case (masked) -> 4
    return A


def best_in_window(A, P, lo, hi):
    """max log-odds of PWM P (len x 4) over start depths lo..hi; A = (n, L), depth d -> col L-d."""
    L = A.shape[1]
    ln = len(P)
    Pext = np.concatenate([P, np.full((ln, 1), -99.0)], axis=1)
    best = np.full(len(A), -99.0)
    cols = np.arange(ln)[None, :]
    for d in range(max(lo, ln), min(hi, L) + 1):
        j = L - d
        best = np.maximum(best, Pext[cols, A[:, j:j + ln]].sum(1))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--label", default="C1_transcribed")
    ap.add_argument("--cohort", default="C1")
    ap.add_argument("--exclude-class", default="pseudogene")
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--oe-specific", type=float, default=2.5)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--tag", default="_v2", help="motif run tag (files hpsad_motifs_<label>_zones<tag>.tsv)")
    ap.add_argument("--sources", default="", help="cohort = only these anchor sources, e.g. EPDnew,CAGE; must match "
                                                    "the motif run's cohort, because thresholds are calibrated on it")
    a = ap.parse_args()
    root = Path(a.root)
    mdir, out = root / "datasets" / "motifs", root / "datasets" / "coverage"
    stem = a.label if a.tag == "_v2" else a.label + a.tag          # v2 file names kept as they were
    out.mkdir(parents=True, exist_ok=True)
    csv.field_size_limit(10 ** 8)

    motifs = [r for r in csv.DictReader(open(mdir / f"hpsad_motifs_{a.label}_zones{a.tag}.tsv", encoding="utf-8"), delimiter="\t")
              if r["pass_fdr"] == "True"]
    el_of, el_name = {}, {}
    for r in csv.DictReader(open(mdir / f"hpsad_elements_{a.label}{a.tag}.tsv", encoding="utf-8"), delimiter="\t"):
        el_name[r["element"]] = r["best_jaspar"].split("\t")[-1] + " " + r["representative"]
        for m in r["members"].split(", "):
            el_of[m.split()[0]] = r["element"]

    rows, seqs = [], []
    with gzip.open(root / "datasets" / "cohorts" / "hpsad_segments_v7.tsv.gz", "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cohort"] != a.cohort or r["has_N"] != "0" or int(r["seg_len"]) < 50:
                continue
            if a.exclude_class and r["gene_class"] == a.exclude_class:
                continue
            if a.primary_only and r["anchor_type"] != "primary":
                continue
            if a.sources and r["anchor_source"] not in a.sources.split(","):
                continue
            rows.append(r)
            seqs.append(r["sequence_softmasked"])
    n = len(seqs)
    print(f"{a.label}: {n:,} promoters, {len(motifs)} motifs passing FDR")

    use, el_oe = [], {}
    for m in motifs:
        mid = f"{m['zone']}:M{int(m['rank']):02d}"
        oe, inst = float(m["ratio"]), int(m["instances"])
        if oe > 20 and inst < 60:
            print(f"  skip copy artefact {mid} {m['consensus']} (O/E {oe}, {inst} carriers)")
            continue
        e = el_of.get(mid, mid)
        el_oe[e] = max(el_oe.get(e, 0), oe)
        use.append((mid, e, m))
    span = max(int(m["start_depth"]) + int(m["window"]) // 2 + 3 + len(m["pwm_ACGT"].split(";")) for _, _, m in use)
    L = int(np.ceil(span / 50) * 50)
    rng = np.random.default_rng(a.seed)

    def shuf(s):
        s = s[-L:]
        return "".join(M.dinuc_shuffle(s[i:i + 50], rng) if len(s[i:i + 50]) > 2 else s[i:i + 50]
                       for i in range(0, len(s), 50))
    R, S = encode(seqs, L), encode([shuf(s) for s in seqs], L)
    elements = sorted(set(e for _, e, _ in use), key=lambda e: -el_oe[e])
    real = {e: np.zeros(n, bool) for e in elements}
    shuf_ = {e: np.zeros(n, bool) for e in elements}
    spec_r, spec_s = np.zeros(n, bool), np.zeros(n, bool)      # specificity judged per MOTIF
    spec_el = {}
    for mid, e, m in use:
        P = np.log2((np.array([[float(x) for x in c.split(",")] for c in m["pwm_ACGT"].split(";")]) + 0.01) / 1.04 / 0.25)
        d0, w, inst = int(m["start_depth"]), int(m["window"]), int(m["instances"])
        lo, hi = d0 - w // 2 - 3, d0 + w // 2 + 3
        br, bs = best_in_window(R, P, lo, hi), best_in_window(S, P, lo, hi)
        thr = np.sort(br)[::-1][min(inst, n) - 1]
        real[e] |= br >= thr
        shuf_[e] |= bs >= thr
        if float(m["ratio"]) >= a.oe_specific:
            spec_r |= br >= thr
            spec_s |= bs >= thr
            spec_el.setdefault(e, []).append((mid, br >= thr, bs >= thr))

    # an element "explains" a promoter only through its specific motifs (O/E >= --oe-specific)
    specific = [e for e in elements if e in spec_el]
    cls = np.array([cpg_class(s) for s in seqs])
    Xr = np.array([np.any([r for _, r, _ in spec_el[e]], axis=0) for e in specific]).T if specific else np.zeros((n, 0), bool)
    Xs = np.array([np.any([q for _, _, q in spec_el[e]], axis=0) for e in specific]).T if specific else np.zeros((n, 0), bool)
    expl_r, expl_s = spec_r, spec_s

    summary = {"label": a.label, "promoters": n, "oe_specific": a.oe_specific, "scan_length": L,
               "elements": {e: {"name": el_name.get(e, e), "max_oe": el_oe[e], "specific": e in specific, "specific_real": int(np.any([r for _, r, _ in spec_el[e]], axis=0).sum()) if e in spec_el else 0,
                                "specific_shuffled": int(np.any([q for _, _, q in spec_el[e]], axis=0).sum()) if e in spec_el else 0,
                                "real": int(real[e].sum()), "shuffled": int(shuf_[e].sum())} for e in elements},
               "explained_real": float(expl_r.mean()), "explained_shuffled": float(expl_s.mean()),
               "coverage": float(expl_r.mean() - expl_s.mean()), "by_cpg": {},
               "n_specific_real": np.bincount(Xr.sum(1), minlength=6)[:6].tolist(),
               "n_specific_shuffled": np.bincount(Xs.sum(1), minlength=6)[:6].tolist(),
               "top_combinations": Counter(" + ".join(el_name.get(e, e).split()[0] for e, b in zip(specific, row) if b) or "(none)"
                                           for row in Xr).most_common(20)}
    for c in ("cpg", "noncpg"):
        k = cls == c
        summary["by_cpg"][c] = {"promoters": int(k.sum()), "explained_real": float(expl_r[k].mean()),
                                "explained_shuffled": float(expl_s[k].mean()),
                                "coverage": float(expl_r[k].mean() - expl_s[k].mean())}
        ids = [rows[i]["anchor_id"] for i in np.where(k & ~expl_r)[0]]
        (out / f"ids_{stem}_unexplained_{c}.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")
        summary["by_cpg"][c]["unexplained_written"] = len(ids)
    src = np.array([r["anchor_source"] for r in rows])
    summary["by_anchor_source"] = {}
    for c in sorted(set(src)):
        k = src == c
        summary["by_anchor_source"][c] = {"promoters": int(k.sum()), "explained_real": float(expl_r[k].mean()),
                                          "explained_shuffled": float(expl_s[k].mean()),
                                          "coverage": float(expl_r[k].mean() - expl_s[k].mean())}
    json.dump(summary, open(out / f"coverage_{stem}_summary.json", "w"), indent=1)
    with open(out / f"coverage_{stem}.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["anchor_id", "gene_id", "symbol", "gene_class", "cpg_class", "explained", "n_specific", "elements"])
        for i, r in enumerate(rows):
            present = [e for e in elements if real[e][i]]
            w.writerow([r["anchor_id"], r["gene_id"], r["symbol"], r["gene_class"], cls[i], int(expl_r[i]),
                        int(Xr[i].sum()), ";".join(present)])

    print(f"\n{'element':<8}{'name':<34}{'O/E':>6}{'real':>8}{'shuffled':>10}{'spec.real':>11}{'spec.shuf':>11}")
    for e in elements:
        d = summary["elements"][e]
        print(f"{e:<8}{d['name'][:33]:<34}{d['max_oe']:>6.1f}{d['real']:>8}{d['shuffled']:>10}{d['specific_real']:>11}{d['specific_shuffled']:>11}")
    print(f"\nexplained (>= 1 specific element): real {expl_r.mean():.1%}, shuffled {expl_s.mean():.1%}, "
          f"coverage {summary['coverage']:.1%}")
    for c, d in summary["by_cpg"].items():
        print(f"  {c:<7} {d['promoters']:>6} promoters: real {d['explained_real']:.1%}, shuffled {d['explained_shuffled']:.1%}, "
              f"coverage {d['coverage']:.1%}; unexplained written {d['unexplained_written']:,}")
    for c, d in summary["by_anchor_source"].items():
        print(f"  anchor {c:<10} {d['promoters']:>6} promoters: real {d['explained_real']:.1%}, "
              f"shuffled {d['explained_shuffled']:.1%}, coverage {d['coverage']:.1%}")
    print(f"\nwritten {out}")


if __name__ == "__main__":
    main()
