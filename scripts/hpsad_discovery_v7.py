#!/usr/bin/env python3
"""
hpsad_discovery_v7.py  -  HPSAD step 4: positional k-mer discovery + test bench

STATISTIC (per k-mer, compared only with itself)
  For each k-mer s and upstream position d (k-mer start at -d), the weighted count
  c_s(d) = sum over promoters of weight (1 / family size) where s starts at -d.
  Occurrences touching a repeat-masked base are not counted.
  A positional peak is a window of width w centred on d whose count O is high
  compared with the same k-mer's density in its own flanks (150 nt each side,
  10 nt gap):  E = flank density x w.   Score = -log10 P(Poisson(E) >= O).
  Composition (GC/CpG gradients, base content) cancels because s is compared
  with itself a few hundred nt away; only a localised excess scores.
  Widths w = 1, 5, 11, 21 nt let sharp and fuzzy-position elements both score.
  Each k-mer keeps its best peak.

NULL AND FDR
  The same cohort after a dinucleotide-preserving shuffle inside 100-nt blocks
  (keeps local composition and repeat-mask runs, destroys positioned elements).
  FDR(t) = #null k-mers scoring >= t / #real k-mers scoring >= t, per k.

TEST BENCH (--bench)
  1. Planted motifs: invented 10-nt motifs, each instance drawn with 85% fidelity
     per base, planted into a known fraction of C1 promoters at a known position
     with Gaussian jitter. Recall = planted motif recovered at FDR <= 5%, at the
     right position.
  2. Known core-promoter elements on the real data (positive controls).

Usage (from the project root):
    python scripts\\hpsad_discovery_v7.py --bench          test bench
    python scripts\\hpsad_discovery_v7.py                  discovery on C1
"""

import argparse
import csv
import gzip
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.special import gammainc

L = 2000
WIDTHS = (1, 5, 11, 21)
FLANK, GAP = 150, 10
MIN_O = 8.0
MIN_RATIO = 1.5
KS = (6, 7, 8)
FDR_MAX = 0.05
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
KNOWN = [  # name, IUPAC consensus, expected start window (upstream positions)
    ("TATA box", "TATAWAWR", (22, 40)),
    ("BRE upstream", "SSRCGCC", (30, 45)),
    ("CCAAT box", "CCAAT", (50, 130)),
    ("GC box (SP1)", "GGGCGG", (30, 150)),
]


# ── data ─────────────────────────────────────────────────────────────────────
def load_cohort(root, cohort="C1"):
    csv.field_size_limit(10 ** 7)
    w = {}
    with open(root / "datasets" / "families" / "hpsad_families_v7.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            w[r["anchor_id"]] = float(r["weight"])
    ids, seqs, wts, cls = [], [], [], []
    with gzip.open(root / "datasets" / "cohorts" / "hpsad_segments_v7.tsv.gz", "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cohort"] != cohort or r["has_N"] != "0" or int(r["seg_len"]) != L:
                continue
            ids.append(r["anchor_id"])
            seqs.append(r["sequence_softmasked"])
            wts.append(w[r["anchor_id"]])
            cls.append(r["gene_class"])
    return ids, seqs, np.array(wts, dtype=np.float64), cls


def encode(seqs, mask=True):
    """(n, L) uint8: 0-3 = ACGT, 4 = masked/other."""
    lut = np.full(256, 4, dtype=np.uint8)
    for ch, v in BASE.items():
        lut[ord(ch)] = v
        if not mask:
            lut[ord(ch.lower())] = v
    buf = np.frombuffer("".join(seqs).encode(), dtype=np.uint8)
    return lut[buf].reshape(len(seqs), L)


# ── counting ─────────────────────────────────────────────────────────────────
def profiles(S, w, k):
    """weighted counts (4^k, L): column d-1 = k-mer starting at upstream position d."""
    n = S.shape[0]
    J = L - k + 1                       # start indices j = 0..L-k ; d = L - j
    code = np.zeros((n, J), dtype=np.int64)
    bad = np.zeros((n, J), dtype=bool)
    for t in range(k):
        col = S[:, t:t + J]
        bad |= col == 4
        code = (code << 2) | (col & 3)
    d_idx = (L - np.arange(J)) - 1      # column index d-1
    P = np.zeros((4 ** k) * L, dtype=np.float64)
    ww = np.broadcast_to(w[:, None], (n, J))
    flat = code * L + d_idx[None, :]
    P += np.bincount(flat[~bad], weights=ww[~bad], minlength=(4 ** k) * L)
    return P.reshape(4 ** k, L).astype(np.float32)


# ── scan ─────────────────────────────────────────────────────────────────────
def window_sum(cs, lo, hi):
    """sum over columns [lo, hi) using cumulative sums cs (rows, L+1); lo/hi arrays."""
    lo = np.clip(lo, 0, L)
    hi = np.clip(hi, 0, L)
    return cs[:, hi] - cs[:, lo], np.maximum(hi - lo, 0)


def scan(P, k, chunk=4096):
    """best peak per k-mer: arrays score, d (centre, upstream position), width, O, E."""
    nk = P.shape[0]
    best = {key: np.zeros(nk, dtype=np.float64) for key in ("score", "d", "w", "O", "E")}
    centres = np.arange(L)                                   # column index = d-1
    valid = centres >= k - 1                                  # k-mer needs d >= k
    for a in range(0, nk, chunk):
        cs = np.zeros((min(chunk, nk - a), L + 1), dtype=np.float64)
        cs[:, 1:] = np.cumsum(P[a:a + chunk], axis=1)
        bs = np.zeros(cs.shape[0])
        bd, bw, bo, be = (np.zeros(cs.shape[0]) for _ in range(4))
        for w in WIDTHS:
            h = w // 2
            O, _ = window_sum(cs, centres - h, centres + h + 1)
            sl, nl = window_sum(cs, centres - h - GAP - FLANK, centres - h - GAP)
            sr, nr = window_sum(cs, centres + h + 1 + GAP, centres + h + 1 + GAP + FLANK)
            dens = (sl + sr + 0.5) / np.maximum(nl + nr, 1)
            E = dens * w
            ok = (O >= MIN_O) & (O >= MIN_RATIO * E) & valid[None, :] & ((nl + nr) >= FLANK)
            with np.errstate(divide="ignore"):
                sc = np.where(ok, -np.log10(np.maximum(gammainc(np.maximum(O, 1e-9), E), 1e-300)), 0.0)
            i = sc.argmax(axis=1)
            s = sc[np.arange(sc.shape[0]), i]
            up = s > bs
            bs[up], bd[up], bw[up] = s[up], i[up] + 1, w
            bo[up] = O[np.arange(O.shape[0]), i][up]
            be[up] = E[np.arange(E.shape[0]), i][up]
        sl_ = slice(a, a + cs.shape[0])
        best["score"][sl_], best["d"][sl_], best["w"][sl_], best["O"][sl_], best["E"][sl_] = bs, bd, bw, bo, be
    return best


def decode(c, k):
    return "".join("ACGT"[(c >> (2 * (k - 1 - t))) & 3] for t in range(k))


def fdr_threshold(real, null):
    """smallest score t with FDR(t) <= FDR_MAX (inf if none)."""
    r = np.sort(real[real > 0])[::-1]
    nl = np.sort(null[null > 0])[::-1]
    best_t = np.inf
    for t in r:
        n_r = np.searchsorted(-r, -t, side="right")
        n_n = np.searchsorted(-nl, -t, side="right")
        if n_n / n_r <= FDR_MAX:
            best_t = t
    return best_t


# ── shuffle (dinucleotide-preserving, Altschul-Erickson) ────────────────────
def dinuc_shuffle(s, rng):
    if len(s) < 3:
        return s
    edges = defaultdict(list)
    for a, b in zip(s, s[1:]):
        edges[a].append(b)
    last = s[-1]
    verts = [v for v in edges if v != last]
    ok = False
    for _ in range(50):
        chosen = {v: int(rng.integers(len(edges[v]))) for v in verts}
        ok = True
        for v in verts:
            u, seen = v, set()
            while u != last:
                if u in seen or u not in chosen:
                    ok = False
                    break
                seen.add(u)
                u = edges[u][chosen[u]]
            if not ok:
                break
        if ok:
            break
    if not ok:      # always valid: the edge leaving each symbol's last occurrence
        chosen = {v: len(edges[v]) - 1 for v in verts}
    lists = {}
    for v, lst in edges.items():
        lst = list(lst)
        if v in chosen:
            tail = lst.pop(chosen[v])
            rng.shuffle(lst)
            lst.append(tail)
        else:
            rng.shuffle(lst)
        lists[v] = lst
    ptr = defaultdict(int)
    out = [s[0]]
    for _ in range(len(s) - 1):
        u = out[-1]
        out.append(lists[u][ptr[u]])
        ptr[u] += 1
    return "".join(out)


def shuffle_cohort(seqs, rng, block=100):
    return ["".join(dinuc_shuffle(s[i:i + block], rng) for i in range(0, len(s), block)) for s in seqs]


# ── run all k ────────────────────────────────────────────────────────────────
def run(seqs, w, null_seqs, log=print):
    out = {}
    S, Sn = encode(seqs), encode(null_seqs)
    for k in KS:
        t0 = time.time()
        real = scan(profiles(S, w, k), k)
        null = scan(profiles(Sn, w, k), k)
        t = fdr_threshold(real["score"], null["score"])
        sig = np.flatnonzero(real["score"] >= t)
        log(f"  k={k}: {len(sig)} k-mers at FDR <= {FDR_MAX} (score >= {t:.1f}); {time.time() - t0:.0f} s")
        out[k] = {"real": real, "null": null, "threshold": t, "sig": sig}
    return out


def hits_table(res):
    rows = []
    for k, r in res.items():
        for c in r["sig"]:
            rows.append({"k": k, "kmer": decode(int(c), k), "d": int(r["real"]["d"][c]),
                         "width": int(r["real"]["w"][c]), "O": round(float(r["real"]["O"][c]), 1),
                         "E": round(float(r["real"]["E"][c]), 1),
                         "ratio": round(float(r["real"]["O"][c] / max(r["real"]["E"][c], 1e-9)), 2),
                         "score": round(float(r["real"]["score"][c]), 1)})
    rows.sort(key=lambda x: -x["score"])
    return rows


# ── test bench ───────────────────────────────────────────────────────────────
def random_motif(rng, length=10):
    while True:
        m = "".join(rng.choice(list("ACGT"), length))
        gc = (m.count("G") + m.count("C")) / length
        if 0.3 <= gc <= 0.7 and max(m.count(b) for b in "ACGT") <= 5 and "CG" not in m:
            return m


def plant(seqs, rng, fractions=(0.005, 0.01, 0.02, 0.05), positions=(30, 90, 250, 700, 1500),
          jitters=(0, 5, 15), fidelity=0.85):
    seqs = [list(s) for s in seqs]
    n = len(seqs)
    scen = []
    for p in positions:
        order = rng.permutation(n)          # promoters for this position, disjoint across scenarios
        used = 0
        for f in fractions:
            for j in jitters:
                m = random_motif(rng)
                cnt = int(round(f * n))
                idx = order[used:used + cnt]
                used += cnt
                for g in idx:
                    off = int(round(rng.normal(0, j))) if j else 0
                    off = max(-2 * j, min(2 * j, off))
                    start = L - (p + off)            # motif starts at upstream position p+off
                    inst = [b if rng.random() < fidelity else rng.choice([x for x in "ACGT" if x != b]) for b in m]
                    if 0 <= start and start + len(m) <= L:
                        seqs[g][start:start + len(m)] = inst
                scen.append({"motif": m, "position": p, "fraction": f, "jitter": j, "n_planted": cnt})
    return ["".join(s) for s in seqs], scen


def hamming_le1(a, b):
    return sum(x != y for x, y in zip(a, b)) <= 1


def evaluate(hits, scen):
    for s in scen:
        m, p, j = s["motif"], s["position"], s["jitter"]
        found = []
        for h in hits:
            k = h["k"]
            for o in range(len(m) - k + 1):
                if hamming_le1(h["kmer"], m[o:o + k]) and abs(h["d"] - (p - o)) <= 2 * j + h["width"] // 2 + 3:
                    found.append(h)
                    break
        s["recovered"] = bool(found)
        s["best_hit"] = max(found, key=lambda x: x["score"]) if found else None
    return scen


def iupac_match(kmer, cons):
    """kmer matches a window of the IUPAC consensus (kmer may be shorter or longer)."""
    def fits(a, c):
        return all(x in IUPAC[y] for x, y in zip(a, c))
    if len(kmer) <= len(cons):
        return any(fits(kmer, cons[o:o + len(kmer)]) for o in range(len(cons) - len(kmer) + 1))
    return any(fits(kmer[o:o + len(cons)], cons) for o in range(len(kmer) - len(cons) + 1))


def known_check(hits):
    out = []
    for name, cons, (lo, hi) in KNOWN:
        m = [h for h in hits if iupac_match(h["kmer"], cons) and lo <= h["d"] <= hi]
        out.append({"element": name, "consensus": cons, "window": [lo, hi], "n_hits": len(m),
                    "top": sorted(m, key=lambda x: -x["score"])[:3]})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--fractions", default="0.005,0.01,0.02,0.05")
    ap.add_argument("--jitters", default="0,5,15")
    ap.add_argument("--fidelity", type=float, default=0.85)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    root = Path(a.root)
    out = root / "datasets" / "discovery"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    ids, seqs, w, cls = load_cohort(root)
    print(f"C1: {len(seqs):,} full-length promoters, {w.sum():,.0f} independent")
    t0 = time.time()
    null_seqs = shuffle_cohort(seqs, rng)
    print(f"Shuffled null built ({time.time() - t0:.0f} s)")

    print("Real data")
    res = run(seqs, w, null_seqs)
    hits = hits_table(res)
    known = known_check(hits)
    with open(out / "hpsad_discovery_C1_hits.tsv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(hits[0]) if hits else ["k"], delimiter="\t")
        wr.writeheader()
        wr.writerows(hits)
    summary = {"cohort": "C1", "promoters": len(seqs), "independent": round(float(w.sum()), 1),
               "parameters": {"widths": WIDTHS, "flank": FLANK, "gap": GAP, "min_O": MIN_O,
                              "min_ratio": MIN_RATIO, "ks": KS, "fdr": FDR_MAX, "seed": a.seed},
               "per_k": {k: {"threshold": float(r["threshold"]), "significant": int(len(r["sig"]))}
                         for k, r in res.items()},
               "known_elements": known}

    if a.bench:
        print("Test bench: planted motifs")
        fr = tuple(float(x) for x in a.fractions.split(","))
        jt = tuple(int(x) for x in a.jitters.split(","))
        planted, scen = plant(seqs, rng, fractions=fr, jitters=jt, fidelity=a.fidelity)
        res_p = run(planted, w, null_seqs)
        hits_p = hits_table(res_p)
        scen = evaluate(hits_p, scen)
        # hits in the planted run explained neither by a planted motif nor by the real data
        real_set = {(h["kmer"]) for h in hits}
        expl = set()
        for sc in scen:
            m = sc["motif"]
            for h in hits_p:
                if any(hamming_le1(h["kmer"], m[o:o + h["k"]]) for o in range(len(m) - h["k"] + 1)):
                    expl.add(h["kmer"])
        unexplained = [h for h in hits_p if h["kmer"] not in real_set and h["kmer"] not in expl]
        summary["bench"] = scen
        summary["bench_settings"] = {"fractions": fr, "jitters": jt, "fidelity": a.fidelity}
        summary["bench_hits"] = len(hits_p)
        summary["bench_unexplained_hits"] = len(unexplained)
        summary["bench_unexplained_examples"] = unexplained[:20]
        print(f"Planted run: {len(hits_p)} hits; {len(unexplained)} explained neither by a planted motif nor by the real data")
        rec = defaultdict(list)
        for s in scen:
            rec[(s["fraction"], s["jitter"])].append(s["recovered"])
        print("Recall by fraction x jitter (over 5 positions):")
        for key in sorted(rec):
            print(f"  fraction {key[0]:.1%} ({int(round(key[0] * len(seqs)))} promoters)  jitter {key[1]:>2} nt : {sum(rec[key])}/{len(rec[key])}")

    with open(out / (f"hpsad_bench{a.tag}_summary.json" if a.bench else "hpsad_discovery_C1_summary.json"),
              "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, default=float)
    print("Known elements:")
    for kn in known:
        top = kn["top"][0] if kn["top"] else None
        print(f"  {kn['element']:<14} {kn['n_hits']:>3} hits" +
              (f"; best {top['kmer']} at -{top['d']} (w={top['width']}, O/E={top['ratio']}, score {top['score']})" if top else ""))
    print(f"Top 25 hits:")
    for h in hits[:25]:
        print(f"  {h['kmer']:<9} -{h['d']:<5} w={h['width']:<3} O={h['O']:<7} E={h['E']:<6} x{h['ratio']:<5} score {h['score']}")


if __name__ == "__main__":
    main()
