#!/usr/bin/env python3
"""
hpsad_syntax_v7.py  -  layout rules between elements: strand bias, and who fills the TATA slot

Uses the same elements, representative PWMs and calibrated thresholds as hpsad_profile_v7.py.

1. Strand bias. In each element's main summit window (11 nt) and over the core (-200..-1):
   sense share of the EXCESS matches, (Os - Es) / ((Os - Es) + (Oa - Ea)), next to the sense share
   expected from the shuffles, Es / (Es + Ea). Palindromic elements (PWM similar to its own reverse
   complement) match both strands at one site and are flagged.
2. TATA slot. Each promoter is called TATA+ when either TATA motif (sense strand, calibrated threshold)
   starts at -37..-25, TATA- otherwise. Every element is profiled separately in the two groups, each
   against the shuffled versions of the same promoters (composition controlled within each group).
   Two slots are compared:  slot A = starts -58..-46,  slot B (the TATA position) = starts -37..-18.
   Promoters are also split by CpG class (-500..-1: GC >= 0.55 and CpG o/e >= 0.6).

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_syntax_v7.py --label C1_evidence --tags _v5s1,_v5s2,_v5s3
    python scripts\\hpsad_syntax_v7.py --write-ids      (TATA group id lists for discovery within groups)
Outputs (datasets\\profiles):
    syntax_strand_<label><tag1>.tsv      strand bias per element
    syntax_tata_slot_<label><tag1>.tsv   O, E, O/E per element x group x slot
    syntax_tata_profile_<label><tag1>.tsv  O and E per position for TATA+ and TATA- (for the figure)
"""
import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402
import hpsad_profile_v7 as Pf                                 # noqa: E402

L = M.L
SLOTS = {"A (-58..-46)": (46, 58), "B, TATA position (-37..-18)": (18, 37)}
TATA_CALL = (25, 37)


def cpg_class(seq):
    u = seq[-500:].upper().replace("N", "")
    n = max(len(u), 1)
    c, g, cg = u.count("C"), u.count("G"), u.count("CG")
    oe = cg * n / (c * g) if c and g else 0.0
    return "cpg" if (c + g) / n >= 0.55 and oe >= 0.6 else "noncpg"


def hits(S, P, thr):
    """(either-strand, sense, antisense) bool matrices; column j <-> start at upstream position L - j."""
    hs = Pf.scan_hits(S, P, thr)
    ha = Pf.scan_hits(S, Pf.rc(P), thr)
    return hs | ha, hs, ha


def by_depth(v):
    return v[::-1]                        # index i <-> upstream position m + i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--label", default="C1_evidence")
    ap.add_argument("--tags", default="_v5s1,_v5s2,_v5s3")
    ap.add_argument("--cohort", default="C1")
    ap.add_argument("--gene-class", default="!pseudogene")
    ap.add_argument("--sources", default="EPDnew,CAGE")
    ap.add_argument("--min-reps", type=int, default=2)
    ap.add_argument("--nulls", type=int, default=10)
    ap.add_argument("--families", default="hpsad_families_v7_down100.tsv", help="family weights file in datasets/families")
    ap.add_argument("--suffix", default="", help="appended to output names (e.g. _fam2)")
    ap.add_argument("--write-ids", action="store_true",
                    help="write anchor-id lists of the TATA groups (for discovery within groups) and stop")
    a = ap.parse_args()
    root = Path(a.root)
    mdir, out = root / "datasets" / "motifs", root / "datasets" / "profiles"
    out.mkdir(parents=True, exist_ok=True)
    tags = a.tags.split(",")
    ids, seqs, w, cls = M.load_cohort(root, a.cohort, a.gene_class, False,
                                      sources=a.sources.split(",") if a.sources else None, families=a.families)
    S = M.encode(seqs)
    shuf = [M.encode(gzip.open(mdir / "nulls" / f"null50_{a.label}{tags[0]}_rep{r}.txt.gz", "rt").read().split("\n"))
            for r in range(a.nulls)]
    cpg = np.array([cpg_class(s) for s in seqs])
    els = Pf.build_elements(mdir, a.label, tags, a.min_reps)
    for el in els:
        rep = el["rep"]
        el["P"] = np.log2((rep["pwm"] + 0.01) / 1.04 / 0.25)
        el["thr"] = Pf.calibrate(S, el["P"], rep["d"], rep["w"], rep["n"])
        el["pal"] = Pf.similar_any_pos(dict(rep), dict(rep, pwm=Pf.rc(rep["pwm"])))

    # ── TATA call (sense strand, starts -37..-25), from the two TATA motifs
    tata = np.zeros(len(seqs), bool)
    tata_els = [el for el in els if "AAA" in el["rep"]["cons"] and el["rep"]["d"] <= 40]     # WWWAAARN, YATAAAN
    for el in tata_els:
        hs = Pf.scan_hits(S, el["P"], el["thr"])
        m = el["P"].shape[0]
        cols = [L - d for d in range(TATA_CALL[0], TATA_CALL[1] + 1) if L - d <= L - m]
        tata |= hs[:, cols].any(axis=1)
    groups = {"TATA+": tata, "TATA-": ~tata, "TATA- CpG": ~tata & (cpg == "cpg"), "TATA- non-CpG": ~tata & (cpg == "noncpg"),
              "TATA+ non-CpG": tata & (cpg == "noncpg"), "all": np.ones(len(seqs), bool)}
    print(f"{a.label}: {len(seqs):,} promoters; TATA+ {tata.sum():,} ({tata.mean():.1%}; "
          f"CpG {(tata & (cpg == 'cpg')).sum():,}, non-CpG {(tata & (cpg == 'noncpg')).sum():,}); "
          f"CpG-island-like {(cpg == 'cpg').sum():,}")
    print("TATA motifs used for the call: " + ", ".join(el["rep"]["cons"] for el in tata_els))
    if a.write_ids:
        for g, fn in (("TATA+", "tata_pos"), ("TATA- CpG", "tata_neg_cpg"), ("TATA- non-CpG", "tata_neg_noncpg")):
            path = out / f"ids_{a.label}_{fn}.txt"
            path.write_text("\n".join(i for i, k in zip(ids, groups[g]) if k) + "\n", encoding="utf-8")
            print(f"  {g:<14} {int(groups[g].sum()):>6,} ids -> {path}")
        return

    strand_rows, slot_rows, prof_rows = [], [], []
    for el in els:
        rep, P = el["rep"], el["P"]
        m = P.shape[0]
        depths = np.arange(m, L + 1)
        ha_, hs_, han_ = hits(S, P, el["thr"])
        sh = [hits(Z, P, el["thr"]) for Z in shuf]
        # strand bias (all promoters)
        Os, Oa = by_depth(w @ hs_), by_depth(w @ han_)
        Es = by_depth(np.mean([w @ x[1] for x in sh], axis=0))
        Ea = by_depth(np.mean([w @ x[2] for x in sh], axis=0))
        Oany = by_depth(w @ ha_)
        Eany = by_depth(np.mean([w @ x[0] for x in sh], axis=0))
        core = depths <= 200
        sm = lambda v: Pf.running_sum(v, 11)
        lr = np.log2((sm(Oany) + 0.5) / (sm(Eany) + 0.5))
        i0 = int(np.argmax(np.where(core, lr, -99)))
        win = (depths >= depths[i0] - 5) & (depths <= depths[i0] + 5)
        for region, sel in (("summit", win), ("core -200..-1", core)):
            xs, xa = Os[sel].sum() - Es[sel].sum(), Oa[sel].sum() - Ea[sel].sum()
            strand_rows.append(dict(element=el["name"], consensus=rep["cons"], palindrome=el["pal"], region=region,
                                    summit=-int(depths[i0]), O_sense=round(Os[sel].sum(), 1), O_antisense=round(Oa[sel].sum(), 1),
                                    E_sense=round(Es[sel].sum(), 1), E_antisense=round(Ea[sel].sum(), 1),
                                    sense_share_observed=round(Os[sel].sum() / max(Os[sel].sum() + Oa[sel].sum(), 1e-9), 2),
                                    sense_share_expected=round(Es[sel].sum() / max(Es[sel].sum() + Ea[sel].sum(), 1e-9), 2),
                                    sense_share_of_excess=round(xs / (xs + xa), 2) if xs + xa > 0 else ""))
        # slots by group
        for g, mask in groups.items():
            wg = w * mask
            Og = by_depth(wg @ ha_)
            Eg = by_depth(np.mean([wg @ x[0] for x in sh], axis=0))
            ng = by_depth(ha_[mask].sum(0).astype(float))
            for sname, (lo, hi) in SLOTS.items():
                sel = (depths >= lo) & (depths <= hi)
                o, e = Og[sel].sum(), Eg[sel].sum()
                slot_rows.append(dict(element=el["name"], consensus=rep["cons"], group=g, promoters_in_group=int(mask.sum()),
                                      slot=sname, O=round(o, 1), E=round(e, 1), ratio=round(o / max(e, 1e-9), 2),
                                      score=round(float(M.poisson_score(np.array([o]), np.array([max(e, 1e-9)]))[0]), 1),
                                      promoters_with_match=int(ng[sel].sum())))
            if g in ("TATA+", "TATA-"):
                for i, d in enumerate(depths):
                    if d <= 200:
                        prof_rows.append(dict(element=el["name"], consensus=rep["cons"], group=g, position=-int(d),
                                              O=round(float(Og[i]), 2), E=round(float(Eg[i]), 3)))
        print(f"  {el['name']:<20} {rep['cons']:<18} summit {-int(depths[i0])}{'  (palindrome)' if el['pal'] else ''}")

    stem = f"{a.label}{tags[0]}{a.suffix}"
    for name, rows in (("syntax_strand", strand_rows), ("syntax_tata_slot", slot_rows), ("syntax_tata_profile", prof_rows)):
        with open(out / f"{name}_{stem}.tsv", "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            wr.writeheader()
            wr.writerows(rows)
    print("\nStrand bias at the main summit (sense share of matches: observed / expected; of the excess):")
    for r in strand_rows:
        if r["region"] == "summit":
            print(f"  {r['element']:<20} {r['consensus']:<18} {r['summit']:>5}  {r['sense_share_observed']:.2f} / "
                  f"{r['sense_share_expected']:.2f}   excess {r['sense_share_of_excess']}{'  palindrome' if r['palindrome'] else ''}")
    print("\nSlots by TATA status (O/E, promoters with a match):")
    gl = ["TATA+", "TATA-", "TATA- CpG", "TATA- non-CpG", "TATA+ non-CpG"]
    print(f"  {'element':<20}{'slot':<10}" + "".join(f"{g:>17}" for g in gl))
    for el in els:
        for sname in SLOTS:
            rr = {r["group"]: r for r in slot_rows if r["consensus"] == el["rep"]["cons"] and r["slot"] == sname}
            print(f"  {el['name'][:19]:<20}{sname[:9]:<10}" + "".join(f"{rr[g]['ratio']:>9.2f} ({rr[g]['promoters_with_match']:>4})" for g in gl))
    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
