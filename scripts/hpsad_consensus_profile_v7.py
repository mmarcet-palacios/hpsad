#!/usr/bin/env python3
"""
hpsad_consensus_profile_v7.py  -  sense-strand profiles of published consensus sequences (Fig 6A: TCT, DPE)

Exact IUPAC match on the sense strand starting at each position -70..+100 of the composite sequences
(+1 = start-site base, no position 0), family-weighted, in the main cohort (C1, no pseudogenes,
EPDnew/CAGE anchors) and in the ten shuffled cohorts of the C1ev_down run (E = their mean).
Usage (from C:\\Promoter\\HPSAD_v7):  python scripts\\hpsad_consensus_profile_v7.py --tag _f
Output: datasets/profiles/downstream_consensus_profiles_C1ev_down<tag>.tsv  (element, consensus, position, O, E)
"""
import argparse, csv, gzip, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_downstream_v7 as Dn, hpsad_motifs_v7 as M      # noqa: E402

CONS = [("TCT (consensus)", "YYCTTTYY"), ("DPE (consensus)", "RGWYVT")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--tag", default="_f")
    ap.add_argument("--nulls", type=int, default=10)
    a = ap.parse_args()
    root = Path(a.root)
    ids, seqs, w, _ = M.load_cohort(root, "C1", "!pseudogene", False, sources=["EPDnew", "CAGE"], segments=Dn.SEG,
                                    families="hpsad_families_v7_down100.tsv")
    w = np.asarray(w)
    S = M.encode(seqs)
    sh = [M.encode(gzip.open(root / "datasets" / "motifs" / "nulls" / f"null50_C1ev_down{a.tag}_rep{r}.txt.gz", "rt").read().split("\n"))
          for r in range(a.nulls)]

    def prof(X, pat):
        return {p: float(w[Dn.iupac_any(X, pat, p, p)].sum()) for p in range(-70, 101)
                if p != 0 and Dn.col(p) + len(pat) <= X.shape[1]}
    out = root / "datasets" / "profiles" / f"downstream_consensus_profiles_C1ev_down{a.tag}.tsv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["element", "consensus", "position", "O", "E"])
        for name, pat in CONS:
            O = prof(S, pat)
            E = {p: 0.0 for p in O}
            for X in sh:
                for p, v in prof(X, pat).items():
                    E[p] += v / len(sh)
            for p in sorted(O):
                wr.writerow([name, pat, p, round(O[p], 2), round(E[p], 3)])
    print(f"written {out}")


if __name__ == "__main__":
    main()
