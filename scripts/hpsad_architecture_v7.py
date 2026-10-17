#!/usr/bin/env python3
"""
hpsad_architecture_v7.py  -  promoter architectures: which elements each promoter carries, and how architectures relate

For every promoter of the main cohort (C1 with an EPDnew or CAGE start site, no pseudogenes), each element of the
vocabulary is called present when a match to its representative PWM (threshold calibrated as in hpsad_profile_v7.py)
or to its published consensus starts inside the element's positional zone, taken from the positional profiles
(Figures 5 and 6). The same calls are made on a shuffled copy of every promoter (dinucleotide shuffle within 50-nt
blocks over -300..+100, random phase), which gives the false-call rate of each element and the number of promoters
expected in each architecture by chance.

An architecture is the set of elements a promoter carries. Network: one node per architecture carried by at least
MIN_NODE family-weighted promoters, size = promoters; edges join architectures that differ by one element.

Outputs (datasets/profiles): architecture_calls_C1ev_v7.tsv (promoter x element, 0/1), architecture_nodes_C1ev_v7.tsv,
architecture_edges_C1ev_v7.tsv, architecture_summary_C1ev_v7.json
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402
import hpsad_profile_v7 as Pf                                 # noqa: E402
from hpsad_syntax_v7 import cpg_class                         # noqa: E402

L, SHIFT = M.L, 100
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
# element: (short name, profile-element name, first start, last start [positions, no 0], strand)
# element: (short name, core words identifying its representative motif, first start, last start, strand)
UP = [("TATA", ["TATAAA", "WWWAAA"], -37, -25, "sense"),
      ("CCAAT", ["CCAAT"], -130, -45, "both"),
      ("GC box", ["GGGCGG", "GGGGCGGGG"], -80, -40, "both"),
      ("ETS", ["CCGGAAG", "GGAAGT"], -60, -15, "both"),
      ("CRE", ["TGACGTCA", "TGACG"], -65, -25, "both"),
      ("NRF1", ["GCGCA.GCGC", "TGCGCA"], -65, -20, "both"),
      ("THAP11", ["ACTACA"], -160, -60, "both"),
      ("ZBTB33", ["TCTCGCGAGA", "CTCGCGAG"], -65, -18, "both"),
      ("AP-1", ["TGA.TCA"], -60, -45, "both")]
CONS = [("Inr", "BBCABW", -3, -3), ("TCT", "YYCTTTYY", -3, -1)]
YY1 = ("YY1", ["AAGATGG"], 1, 25)
MIN_NODE = 25


def col(pos):
    d = SHIFT - pos + 1 if pos > 0 else SHIFT - pos
    return L - d


def iupac_hit(S, pat, lo, hi):
    lut = np.zeros((len(pat), 5), bool)
    for c, ch in enumerate(pat):
        for b in IUPAC[ch]:
            lut[c, "ACGT".index(b)] = True
    hit = np.zeros(S.shape[0], bool)
    for p in range(lo, hi + 1):
        if p == 0:
            continue
        j = col(p)
        ok = np.ones(S.shape[0], bool)
        for c in range(len(pat)):
            ok &= lut[c][S[:, j + c]]
        hit |= ok
    return hit


def pwm_hit(S, P, thr, lo, hi, strand):
    """any match starting at positions lo..hi (composite columns)."""
    cols = [col(p) for p in range(lo, hi + 1) if p != 0]
    hs = Pf.scan_hits(S, P, thr)[:, cols].any(1)
    if strand == "both":
        hs |= Pf.scan_hits(S, Pf.rc(P), thr)[:, cols].any(1)
    return hs


def shuffle_core(seqs, rng, block=50, span=400):
    out = []
    for s in seqs:
        head, core = s[:-span], s[-span:]
        ph = int(rng.integers(block))
        parts = [core[:ph]] + [core[i:i + block] for i in range(ph, span, block)]
        out.append(head + "".join(M.dinuc_shuffle(p, rng) if len(p) > 2 else p for p in parts))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--families", default="hpsad_families_v7_down100.tsv")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--tags", default="_v5s1,_v5s2,_v5s3", help="replicate tags of the C1_evidence runs (elements)")
    ap.add_argument("--down-tag", default="_v7", help="tag of the C1ev_down runs (YY1)")
    ap.add_argument("--min-reps", type=int, default=2)
    a = ap.parse_args()
    root = Path(a.root)
    mdir, out = root / "datasets" / "motifs", root / "datasets" / "profiles"
    src = ["EPDnew", "CAGE"]
    ids_u, seqs_u, w_u, _ = M.load_cohort(root, "C1", "!pseudogene", False, sources=src, families=a.families)
    ids, seqs, w, cls = M.load_cohort(root, "C1", "!pseudogene", False, sources=src, families=a.families,
                                      segments="hpsad_segments_v7_down100.tsv.gz")
    assert ids == ids_u, "upstream and composite cohorts differ"
    Su, S = M.encode(seqs_u), M.encode(seqs)
    rng = np.random.default_rng(a.seed)
    Z = M.encode(shuffle_core(seqs, rng))
    print(f"{len(ids):,} promoters ({w.sum():,.1f} independent)")

    els = Pf.build_elements(mdir, "C1_evidence", a.tags.split(","), a.min_reps)
    calls, callz, names = [], [], []
    taken = set()
    for short, words, lo, hi, strand in UP:
        # elements whose representative consensus carries the core word (first word that matches any element wins)
        reps = []
        for wd in words:
            reps = [e for k, e in enumerate(els) if k not in taken and M.word_strand(e["rep"]["cons"], [wd])]
            if reps:
                break
        taken |= {k for k, e in enumerate(els) if any(e is r for r in reps)}
        print(f"  {short:<7} <- " + (", ".join(e["rep"]["cons"] for e in reps) or "no element"))
        hr, hz = np.zeros(len(ids), bool), np.zeros(len(ids), bool)
        for e in reps:
            r = e["rep"]
            pwm = r["pwm"] if M.word_strand(r["cons"], words) >= 0 else r["pwm"][::-1, ::-1]
            P = np.log2((pwm + 0.01) / 1.04 / 0.25)
            thr = Pf.calibrate(Su, P, r["d"], r["w"], r["n"])
            hr |= pwm_hit(S, P, thr, lo, hi, strand)
            hz |= pwm_hit(Z, P, thr, lo, hi, strand)
        calls.append(hr); callz.append(hz); names.append(short)
    for short, pat, lo, hi in CONS:
        calls.append(iupac_hit(S, pat, lo, hi)); callz.append(iupac_hit(Z, pat, lo, hi)); names.append(short)
    rows = [r for r in csv.DictReader(open(mdir / f"hpsad_motifs_C1ev_down_zones{a.down_tag}.tsv", encoding="utf-8"),
                                      delimiter="\t")
            if r["pass_fdr"] == "True" and M.word_strand(r["consensus"], YY1[1]) == 1 and YY1[2] <= int(r["position"]) <= YY1[3] + 10]
    r = max(rows, key=lambda x: float(x["score"]))
    print(f"  YY1     <- {r['consensus']} ({r['position']})")
    P = np.log2((np.array([[float(x) for x in c.split(",")] for c in r["pwm_ACGT"].split(";")]) + 0.01) / 1.04 / 0.25)
    # threshold: calibrated on the whole cohort that found it (as in hpsad_downstream_v7.py)
    thr = Pf.calibrate(S, P, int(r["start_depth"]), int(r["window"]), int(r["instances"]))
    calls.append(pwm_hit(S, P, thr, YY1[2], YY1[3], "sense")); callz.append(pwm_hit(Z, P, thr, YY1[2], YY1[3], "sense"))
    names.append("YY1")
    C, CZ = np.array(calls).T, np.array(callz).T                # (n, elements)

    per_el = []
    for k, nm in enumerate(names):
        o, e = float(w[C[:, k]].sum()), float(w[CZ[:, k]].sum())
        per_el.append(dict(element=nm, promoters=int(C[:, k].sum()), weighted=round(o, 1), shuffled=round(e, 1),
                           ratio=round(o / max(e, 0.5), 2), true_fraction=round(max(0.0, 1 - e / max(o, 1e-9)), 2)))
        print(f"  {nm:<8} {int(C[:, k].sum()):>6} promoters  weighted {o:8.1f}  shuffled {e:8.1f}  x{o / max(e, .5):5.1f}")

    def arch(row):
        return tuple(n for n, b in zip(names, row) if b)
    A, AZ = Counter(), Counter()
    for i in range(len(ids)):
        A[arch(C[i])] += w[i]
        AZ[arch(CZ[i])] += w[i]
    cpg = np.array([cpg_class(s) == "cpg" for s in seqs_u])
    keys = [arch(C[i]) for i in range(len(ids))]
    nodes = []
    for k, v in sorted(A.items(), key=lambda x: -x[1]):
        if v < MIN_NODE:
            continue
        m = np.array([kk == k for kk in keys])
        nodes.append(dict(arch=" + ".join(k) if k else "(none)", n_elements=len(k), promoters=int(m.sum()), weighted=round(v, 1),
                          shuffled=round(AZ.get(k, 0), 1), ratio=round(v / max(AZ.get(k, 0), 0.5), 2),
                          cpg_fraction=round(float(w[m & cpg].sum() / w[m].sum()), 3),
                          protein_coding_fraction=round(float(w[m & (np.array(cls) == "protein_coding")].sum() / w[m].sum()), 3)))
    keyset = {tuple(n["arch"].split(" + ")) if n["arch"] != "(none)" else (): n["arch"] for n in nodes}
    edges = []
    ks = list(keyset)
    for x in range(len(ks)):
        for y in range(x + 1, len(ks)):
            a_, b_ = set(ks[x]), set(ks[y])
            if len(a_ ^ b_) == 1:
                edges.append(dict(a=keyset[ks[x]], b=keyset[ks[y]], differs=list(a_ ^ b_)[0]))
    tot = float(w.sum())
    summary = dict(promoters=len(ids), independent=round(tot, 1), elements=per_el,
                   none_fraction=round(A[()] / tot, 3), none_fraction_shuffled=round(AZ[()] / tot, 3),
                   cpg_fraction_all=round(float(w[np.array([cpg_class(s) == "cpg" for s in seqs_u])].sum() / tot), 3),
                   in_nodes_fraction=round(sum(n["weighted"] for n in nodes) / tot, 3), n_architectures=len(A), n_nodes=len(nodes))
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "architecture_calls_C1ev_v7.tsv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["anchor_id", "gene_class", "weight"] + names)
        for i in range(len(ids)):
            wr.writerow([ids[i], cls[i], round(float(w[i]), 4)] + [int(b) for b in C[i]])
    for nm, rows_ in (("architecture_nodes_C1ev_v7.tsv", nodes), ("architecture_edges_C1ev_v7.tsv", edges)):
        with open(out / nm, "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows_[0]), delimiter="\t")
            wr.writeheader()
            wr.writerows(rows_)
    json.dump(summary, open(out / "architecture_summary_C1ev_v7.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "elements"}, indent=1))
    for n in nodes[:40]:
        print(f"  {n['weighted']:8.1f}  (shuffled {n['shuffled']:7.1f}, x{n['ratio']:5.2f})  {n['arch']}")


if __name__ == "__main__":
    main()
