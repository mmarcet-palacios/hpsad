#!/usr/bin/env python3
"""
hpsad_downstream_v7.py  -  the start-site region (-60..+100): profiles and a targeted test

Uses the composite sequences of extract_downstream_v7.py (-1900..-1, +1..+100; +1 = start-site base,
no position 0) and the ten shuffled cohorts of the C1ev_down run (dinucleotide shuffles in 50-nt
blocks with random phase). All counts are on the SENSE strand (these elements act in one orientation
relative to the start site) and family-weighted; E = mean of the ten shuffles.

1. Profiles of discovered elements (PWM, threshold calibrated as in hpsad_profile_v7.py on the cohort
   that found it): Initiator, TATA, YY1, start-codon context, 5' splice donor.
2. Targeted test of published core-promoter consensus sequences at their published positions
   (exact IUPAC match, start within the window): promoters matching in the real cohort vs the
   shuffles, overall, by TATA group, and among Initiator promoters (BBCA+1BW at -3).

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_downstream_v7.py
Outputs (datasets\\profiles):
    downstream_profiles_C1ev_down_v7.tsv   element x position (-60..+100): O, E
    downstream_targeted_C1ev_down_v7.tsv   consensus x group: O, E (mean, sd of shuffles), O/E, -log10 P (Poisson), z vs shuffles, promoters
    downstream_profiles_v7.png             the profiles
"""
import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
from scipy.special import gammainc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402
import hpsad_profile_v7 as Pf                                 # noqa: E402
import hpsad_stability_v7 as St                               # noqa: E402

L, SHIFT = M.L, 100
SEG = "hpsad_segments_v7_down100.tsv.gz"
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

# discovered elements: (label, run label, core words, first and last start position). The best-scoring passing
# motif of the run whose consensus contains a core word (either strand; the PWM is turned to the sense strand) and
# whose start lies in the range represents the element.
ELEMENTS = [("TATA box", "C1ev_down", ["TATAAA", "WWWAAA"], -36, -24),
            ("Initiator, TCA+1", "C1ev_down", ["TCA"], -6, 0),
            ("Initiator, CA+1G", "C1ev_down", ["CAG"], -6, 0),
            ("YY1", "C1ev_down", ["AAGATGG"], 1, 30),
            ("start codon (Kozak)", "C1ev_down", ["CCATGG", "ACCATGG"], 5, 80),
            ("5' splice donor", "C1ev_down", ["CAGGTA", "AGGTAAG", "GGTGAG"], 20, 95)]
# published consensus: (name, IUPAC, first allowed start, last allowed start)  [positions, no 0]
TARGETED = [("TATA box (control)", "TATAWAWR", -33, -29),
            ("Initiator BBCA+1BW (control)", "BBCABW", -3, -3),
            ("TCT  YYC+1TTTYY", "YYCTTTYY", -3, -1),
            ("XCPE1 DSGYGGRASM", "DSGYGGRASM", -9, -7),
            ("BREu SSRCGCC", "SSRCGCC", -39, -36),
            ("BREd RTDKKKK", "RTDKKKK", -24, -21),
            ("MTE CSARCSSAACGS", "CSARCSSAACGS", 17, 19),
            ("MTE core SAACGS", "SAACGS", 23, 25),
            ("DPE RGWYV", "RGWYV", 27, 29),
            ("DPE RGWYVT", "RGWYVT", 27, 29),
            ("DCE SI CTTC", "CTTC", 5, 8),
            ("DCE SII CTGT", "CTGT", 15, 18)]


def col(pos):
    """position (no 0) -> column of the composite sequence."""
    d = SHIFT - pos + 1 if pos > 0 else SHIFT - pos
    return L - d


def pos_of_col(j):
    d = L - j
    return -(d - SHIFT) if d > SHIFT else SHIFT + 1 - d


def iupac_any(S, pat, p_lo, p_hi):
    """bool per promoter: exact match of pat (sense) starting at any position p_lo..p_hi."""
    lut = np.zeros((len(pat), 5), bool)
    for c, ch in enumerate(pat):
        for b in IUPAC[ch]:
            lut[c, "ACGT".index(b)] = True
    hit = np.zeros(S.shape[0], bool)
    for p in range(p_lo, p_hi + 1):
        if p == 0:
            continue
        j = col(p)
        ok = np.ones(S.shape[0], bool)
        for c in range(len(pat)):
            ok &= lut[c][S[:, j + c]]
        hit |= ok
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--nulls", type=int, default=10)
    ap.add_argument("--tag", default="_v7", help="discovery run tag of the C1ev_down runs (motif files and nulls)")
    ap.add_argument("--families", default="hpsad_families_v7_down100.tsv", help="family weights file in datasets/families")
    ap.add_argument("--suffix", default="", help="appended to output names (e.g. _fam2)")
    a = ap.parse_args()
    root = Path(a.root)
    mdir, out = root / "datasets" / "motifs", root / "datasets" / "profiles"
    out.mkdir(parents=True, exist_ok=True)
    src = ["EPDnew", "CAGE"]
    ids, seqs, w, _ = M.load_cohort(root, "C1", "!pseudogene", False, sources=src, segments=SEG, families=a.families)
    S = M.encode(seqs)
    shuf = []
    for r in range(a.nulls):
        ns = gzip.open(mdir / "nulls" / f"null50_C1ev_down{a.tag}_rep{r}.txt.gz", "rt").read().split("\n")
        if len(ns) != len(seqs):
            sys.exit(f"null rep{r}: {len(ns)} sequences vs {len(seqs)}")
        shuf.append(M.encode(ns))
    idx = {i: k for k, i in enumerate(ids)}
    groups = {"all": np.ones(len(ids), bool)}
    for g, fn in (("TATA+", "tata_pos"), ("TATA- CpG", "tata_neg_cpg"), ("TATA- non-CpG", "tata_neg_noncpg")):
        m = np.zeros(len(ids), bool)
        for x in open(out / f"ids_C1_evidence_{fn}.txt", encoding="utf-8"):
            if x.strip() in idx:
                m[idx[x.strip()]] = True
        groups[g] = m
    inr = iupac_any(S, "BBCABW", -3, -3)
    groups["Initiator+"] = inr
    groups["Initiator+ TATA-"] = inr & ~groups["TATA+"]
    print(f"{len(ids):,} promoters; " + ", ".join(f"{g} {int(m.sum()):,}" for g, m in groups.items()))

    # ── 1. profiles of discovered elements
    prof_rows, curves = [], []
    used = set()
    for label, run, words, plo, phi in ELEMENTS:
        cand = []
        rows_ = [r for r in csv.DictReader(open(mdir / f"hpsad_motifs_{run}_zones{a.tag}.tsv", encoding="utf-8"), delimiter="\t")
                 if r["pass_fdr"] == "True" and (run, r["zone"], r["rank"]) not in used and plo <= int(r["position"]) <= phi]
        for wd in words:                                # core words in order of preference
            cand = [(float(r["score"]), M.word_strand(r["consensus"], [wd]), r) for r in rows_]
            cand = [c for c in cand if c[1]]
            if cand:
                break
        if not cand:
            print(f"  {label}: no passing motif with {words} at {plo:+d}..{phi:+d} in {run}{a.tag} - skipped")
            continue
        _, st, r = max(cand, key=lambda x: x[0])
        used.add((run, r["zone"], r["rank"]))
        cons = r["consensus"]
        pwm = np.array([[float(x) for x in c.split(",")] for c in r["pwm_ACGT"].split(";")])
        if st < 0:
            pwm = pwm[::-1, ::-1]
        P = np.log2((pwm + 0.01) / 1.04 / 0.25)
        cal = groups["all"] if run == "C1ev_down" else groups["TATA- non-CpG"]
        thr = Pf.calibrate(S[cal], P, int(r["start_depth"]), int(r["window"]), int(r["instances"]))
        hs = Pf.scan_hits(S, P, thr)
        O = w @ hs
        E = np.mean([w @ Pf.scan_hits(Z, P, thr) for Z in shuf], axis=0)
        J = hs.shape[1]
        pos = np.array([pos_of_col(j) for j in range(J)])
        sel = (pos >= -60) & (pos <= 100 - P.shape[0] + 1)
        for p, o, e in zip(pos[sel], O[sel], E[sel]):
            prof_rows.append(dict(element=label, consensus=cons, position=int(p), O=round(float(o), 2), E=round(float(e), 3)))
        top = pos[sel][np.argmax((O[sel] + 0.5) / (E[sel] + 0.5))]
        print(f"  {label:<22} {cons:<15} threshold {thr:5.2f}  best start {top:+d}  "
              f"O/E {((O[sel] + 0.5) / (E[sel] + 0.5)).max():.1f}  matches at best {O[sel].max():.0f}")
        curves.append((label, cons, pos[sel], O[sel], E[sel]))

    # ── 2. targeted consensus test
    t_rows = []
    print(f"\nO/E (promoters, -log10 P Poisson, z vs shuffles)\n{'element':<30}{'window':>12}" + "".join(f"{g:>27}" for g in groups))
    for name, pat, lo, hi in TARGETED:
        hr = iupac_any(S, pat, lo, hi)
        hn = [iupac_any(Z, pat, lo, hi) for Z in shuf]
        line = f"{name:<30}{f'{lo:+d}..{hi:+d}':>12}"
        for g, m in groups.items():
            o = float(w[m & hr].sum())
            e = float(np.mean([w[m & h].sum() for h in hn]))
            sd = float(np.std([w[m & h].sum() for h in hn]))
            sc = float(-np.log10(max(gammainc(o, max(e, 1e-9)), 1e-300))) if o > e else 0.0   # plain Poisson, no gate
            z = (o - e) / max(sd, 1e-9)
            t_rows.append(dict(element=name, consensus=pat, start_from=lo, start_to=hi, group=g, promoters_in_group=int(m.sum()),
                               promoters_matching=int((m & hr).sum()), O=round(o, 1), E=round(e, 1), E_sd=round(sd, 1),
                               ratio=round(o / max(e, 1e-9), 2), score=round(sc, 1), z_vs_shuffles=round(z, 1)))
            line += f"{o / max(e, 1e-9):>8.2f} ({int((m & hr).sum()):>5},{min(sc, 999):>5.1f},{z:>5.1f})"
        print(line)

    for nm, rows in ((f"downstream_profiles_C1ev_down{a.tag}{a.suffix}.tsv", prof_rows), (f"downstream_targeted_C1ev_down{a.tag}{a.suffix}.tsv", t_rows)):
        with open(out / nm, "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            wr.writeheader()
            wr.writerows(rows)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(len(curves), 1, figsize=(8, 1.2 * len(curves) + 0.8), sharex=True)
        for a_, (label, cons, pos, O, E) in zip(np.atleast_1d(ax), curves):
            k = np.ones(3) / 3
            y = (np.convolve(O, k, "same") + 0.5) / (np.convolve(E, k, "same") + 0.5)
            a_.fill_between(pos, 1, y, where=y > 1, color="#3f6fa8", lw=0)
            a_.plot(pos, y, c="#2c5282", lw=0.6)
            a_.axhline(1, c="#888", lw=0.5)
            a_.axvline(0, c="#b5452f", lw=0.8, ls=(0, (3, 2)))
            a_.text(-58, a_.get_ylim()[1] * 0.78, f"{label}  {cons}", fontsize=8, weight="bold")
            a_.tick_params(labelsize=7)
        np.atleast_1d(ax)[-1].set_xlabel("motif start relative to the start site (nt; dashed = between -1 and +1)", fontsize=8)
        plt.tight_layout()
        plt.savefig(out / f"downstream_profiles{a.tag}{a.suffix}.png", dpi=150)
    except ImportError:
        pass
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
