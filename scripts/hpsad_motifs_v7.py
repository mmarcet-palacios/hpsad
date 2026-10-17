#!/usr/bin/env python3
"""
hpsad_motifs_v7.py  -  HPSAD step 4b: iterative positional motif discovery
("find the strongest positioned motif, explain it, mask it, repeat")

For one cohort of full-length promoters (2,000 nt, last base = position -1),
each promoter weighted 1 / family size, repeat-masked bases excluded:

  EXPECTATION. A motif's count O in a window of width w is compared with
            E = w x (its own density in two 150-nt flanks, 10-nt gap)
                  x (window density / flank density of the same motif in a
                     composition reference: the cohort after a dinucleotide-
                     preserving shuffle in 50-nt blocks, random phase).
            The reference keeps the GC/CpG gradient towards the start site but
            no positioned elements, so a pure composition effect cancels.
            score = -log10 P(Poisson(E) >= O).

  1. SEED   Every exact 6-, 7- and 8-mer is scored once for its best positional
            peak (window w in 1..51 nt). Words
            are kept in a queue; at each round the best word is re-scored on
            the current (masked) data and used as a seed only if it is still
            the best (lazy greedy).
  2. REFINE The seed's occurrences around its peak are aligned, extended by
            4 nt on each side and turned into a position weight matrix (PWM).
            The PWM is scanned over the window and its flanks; the match
            threshold, window width and centre that give the strongest
            positional peak are chosen; matches inside the window become the
            new instances. Repeat until the instances stop changing (max 8).
            End columns with < 0.25 bits are trimmed.
  3. MASK   The bases of all instances inside the window are masked, so the
            same signal (and its variants) cannot seed again.
  4. Repeat until no seed scores >= STOP_SCORE or MAX_MOTIFS is reached.

FDR: the identical loop (same composition reference) is run on R further null
cohorts (independent shuffles of the same kind, repeat mask carried along). For motif score u and window class (sharp: w <= 11 nt;
broad: w >= 21 nt):  FDR(u) = mean_r #null motifs >= u / #real motifs >= u.

ZONES: discovery runs separately for the core zone (peak centre at -1..-300) and
the distal zone (-301..-2000), each with its own nulls and FDR, so the chance
ceiling of the core reflects only the core's positions.

Usage (from the project root):
    python scripts\\hpsad_motifs_v7.py --bench        test bench (C1, planted motifs)
    python scripts\\hpsad_motifs_v7.py --cohort C1    discovery
    options: --nulls 3  --workers N (shuffles)  --jobs N (parallel discovery runs,
             ~4 GB RAM each; e.g. --jobs 5 runs 3 nulls + real + planted together)
"""

import argparse
import concurrent.futures as cf
import csv
import gzip
import heapq
import json
import math
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.special import gammainc

# Optional GPU backend for the word scan (set HPSAD_GPU=1, or --gpu in run_cohorts_parallel_v7.py).
# Needs CuPy (pip install "cupy-cuda13x[ctk]"). Same arithmetic in float64; falls back to CPU if missing.
GPU = False
if os.environ.get("HPSAD_GPU", "") == "1":
    try:
        import cupy as cp
        from cupyx.scipy.special import gammainc as cp_gammainc
        cp.zeros(1).sum()
        GPU = True
    except Exception as _e:                                   # no CuPy / no device: stay on CPU
        print(f"HPSAD_GPU=1 but the GPU backend is unavailable ({type(_e).__name__}: {_e}); using CPU")

L = 2000
WIDTHS = (1, 5, 11, 21, 31, 51)
FLANK, GAP = 150, 10
MIN_O = 8.0
MIN_RATIO = 1.5
KS = (6, 7, 8)
POOL_PER_K = 4000         # best words per k kept as seed candidates
GPU_BATCH = int(os.environ.get("HPSAD_GPU_BATCH", "1024"))   # words per GPU call
SCORE_DECIMALS = 9        # window scores are rounded to this many decimals before any comparison
EXT = 4                   # PWM extension on each side of the seed
PAD = 10                  # extra search margin around the window
STOP_SCORE = 10.0
MAX_MOTIFS = 40
ZONES = {"core": (1, 300), "distal": (301, L)}     # peak centre (upstream position) ranges
FDR_MAX = 0.05
TRIM_BITS = 0.25
BASES = "ACGT"
_RNG = np.random.default_rng(0)      # sampling for threshold candidates
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
KNOWN = [("TATA box", "TATAWAW", (20, 42)), ("BRE upstream", "SSRCGCC", (30, 45)),
         ("CCAAT box (NF-Y)", "CCAAT", (30, 140)), ("GC box (SP1)", "GGGCGG", (25, 160)),
         ("ETS", "CCGGAAG", (10, 120)), ("E-box", "CACGTG", (10, 150)),
         ("CRE (CREB/ATF)", "TGACGTCA", (20, 150)), ("NRF1", "GCGCATGCGC", (5, 120)),
         ("YY1", "AAGATGGCG", (1, 40))]


# ── data ─────────────────────────────────────────────────────────────────────
# small-RNA genes: the annotated 5' end is the end of a short product (tRNA, snRNA, snoRNA, miRNA hairpin ...),
# transcribed by Pol III, by snRNA-type promoters or from a host transcript - not a Pol II core promoter.
SMALL_RNA = {"tRNA", "snRNA", "snoRNA", "miRNA", "rRNA", "misc_RNA", "scRNA", "vault_RNA", "Y_RNA",
             "RNase_MRP_RNA", "RNase_P_RNA", "telomerase_RNA"}


def gene_class_of(r):
    """gene class of a segment row; small-RNA genes form their own class 'small_rna'."""
    return "small_rna" if r.get("biotype", "") in SMALL_RNA else r["gene_class"]


def load_cohort(root, cohort="C1", gene_class=None, primary_only=False, min_len=50, ids_file=None, sources=None,
                segments="hpsad_segments_v7.tsv.gz", families="hpsad_families_v7_down100.tsv"):
    """segments of one cohort (optionally one gene class). Segments shorter than L are
    padded upstream with 'n' (masked: never counted, never used as flank evidence).
    gene_class: a class ('protein_coding', 'non_coding', 'pseudogene', 'small_rna'), or '!a,b' = all classes but a and b.
    Note: 'non_coding' no longer includes small-RNA genes (biotypes in SMALL_RNA)."""
    csv.field_size_limit(10 ** 7)
    wmap = {}
    with open(root / "datasets" / "families" / families, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            wmap[r["anchor_id"]] = float(r["weight"])
    keep = None
    if ids_file:                                              # restrict to listed anchor ids
        keep = set(x.strip() for x in open(Path(root) / ids_file, encoding="utf-8") if x.strip())
    ids, seqs, wts, cls = [], [], [], []
    with gzip.open(root / "datasets" / "cohorts" / segments, "rt", encoding="utf-8") as fh:   # or the +D composite
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cohort"] != cohort or r["has_N"] != "0" or int(r["seg_len"]) < min_len:
                continue
            gc = gene_class_of(r)
            if gene_class and gene_class.startswith("!"):          # "!pseudogene" = all but pseudogenes
                if gc in gene_class[1:].split(","):
                    continue
            elif gene_class and gc != gene_class:
                continue
            if primary_only and r["anchor_type"] != "primary":
                continue
            if keep is not None and r["anchor_id"] not in keep:
                continue
            if sources and r.get("anchor_source", "") not in sources:   # e.g. EPDnew, CAGE
                continue
            seq = r["sequence_softmasked"][-L:]
            ids.append(r["anchor_id"])
            seqs.append("n" * (L - len(seq)) + seq)
            wts.append(wmap[r["anchor_id"]])
            cls.append(gc)
    return ids, seqs, np.array(wts), cls


def encode(seqs):
    """(n, L) uint8: 0-3 = ACGT (upper case only), 4 = repeat-masked / other."""
    lut = np.full(256, 4, dtype=np.uint8)
    for i, ch in enumerate(BASES):
        lut[ord(ch)] = i
    return lut[np.frombuffer("".join(seqs).encode(), dtype=np.uint8)].reshape(len(seqs), L)


# ── null cohort ──────────────────────────────────────────────────────────────
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
    if not ok:
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


def shuffle_one(s, rng, block=50):
    phase = int(rng.integers(0, block))
    cuts = [0] + list(range(phase, len(s), block)) + [len(s)]
    cuts = sorted(set(cuts))
    return "".join(dinuc_shuffle(s[a:b], rng) for a, b in zip(cuts, cuts[1:]))


def _shuffle_chunk(args):
    """shuffle the real part of each sequence only; upstream 'n' padding is kept."""
    seqs, seed = args
    rng = np.random.default_rng(seed)
    out = []
    for s in seqs:
        core = s.lstrip("n")
        out.append("n" * (len(s) - len(core)) + (shuffle_one(core, rng) if len(core) >= 3 else core))
    return out


def null_cohort(seqs, seed, path, workers):
    if path.exists():
        with gzip.open(path, "rt") as fh:
            cached = fh.read().split("\n")
        if len(cached) == len(seqs):
            return cached
    n_chunks = max(1, workers) * 4
    size = -(-len(seqs) // n_chunks)
    jobs = [(seqs[i:i + size], seed * 1000 + j) for j, i in enumerate(range(0, len(seqs), size))]
    if workers > 1:
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            parts = list(ex.map(_shuffle_chunk, jobs))
    else:
        parts = [_shuffle_chunk(j) for j in jobs]
    out = [x for p in parts for x in p]
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(out))
    return out


# ── positional scoring of a count profile ────────────────────────────────────
def poisson_score(O, E, xp=np, ginc=gammainc):
    """-log10 P(Poisson(E) >= O); arrays (NumPy, or CuPy with xp=cupy)."""
    sc = xp.zeros(O.shape)
    ok = (O >= MIN_O) & (O >= MIN_RATIO * E)
    if bool(ok.any()):
        sc[ok] = -xp.log10(xp.maximum(ginc(O[ok], E[ok]), 1e-300))
    return sc


def best_peak(cs, lo_d, hi_d, k=1, centre_range=None, csr=None, xp=np, ginc=gammainc, csf=None, ca=None):
    """best window over a count profile.
    cs : cumulative sums (rows, ncol+1) of counts by depth, column 0 = depth lo_d.
    csr: the same for the composition reference cohort (optional). Then
         E = w x own flank density x (reference window density / reference flank density),
         so a pure composition gradient cancels.
    csf: cumulative counts used for the FLANKS (default cs). After masking, flanks are counted on the
         unmasked sequence, so earlier motifs' instances do not deplete the expectation.
    ca : cumulative availability per column (fraction of opportunity left after masking, 0..1);
         the window's expectation is scaled by its summed availability instead of its width.
    Returns per row: score, centre depth, width, O, E."""
    if csf is None:
        csf = cs
    ncol = cs.shape[1] - 1
    allc = np.arange(ncol)
    valid = (allc + lo_d) >= k
    if centre_range is not None:
        valid &= ((allc + lo_d) >= centre_range[0]) & ((allc + lo_d) <= centre_range[1])
    centres = xp.asarray(allc[valid])         # score only centres that can win (same result, less work)
    depth = centres + lo_d
    rows = xp.arange(cs.shape[0])
    if len(centres) == 0:
        return [np.zeros(cs.shape[0]) for _ in range(5)]
    best = [xp.zeros(cs.shape[0]) for _ in range(5)]
    for w in WIDTHS:
        h = w // 2

        def wsum(a, b):
            a = xp.clip(a, 0, ncol)
            b = xp.clip(b, 0, ncol)
            return cs[:, b] - cs[:, a], xp.maximum(b - a, 0)
        def fsum(a, b):
            a = xp.clip(a, 0, ncol)
            b = xp.clip(b, 0, ncol)
            return csf[:, b] - csf[:, a], xp.maximum(b - a, 0)
        O, _ = wsum(centres - h, centres + h + 1)
        sl, nl = fsum(centres - h - GAP - FLANK, centres - h - GAP)
        sr, nr = fsum(centres + h + 1 + GAP, centres + h + 1 + GAP + FLANK)
        nf = xp.maximum(nl + nr, 1)
        if ca is None:
            wa = wn = w
        else:
            a_, b_ = xp.clip(centres - h, 0, ncol), xp.clip(centres + h + 1, 0, ncol)
            wa = ca[:, b_] - ca[:, a_]              # available columns (masking) ...
            wn = xp.maximum(b_ - a_, 1)             # ... out of the columns actually in the window
        E = (sl + sr + 0.5) / nf * wa
        if csr is not None:
            def rsum(a, b):
                a = xp.clip(a, 0, ncol)
                b = xp.clip(b, 0, ncol)
                return csr[:, b] - csr[:, a]
            Or = rsum(centres - h, centres + h + 1)
            Fr = rsum(centres - h - GAP - FLANK, centres - h - GAP) + rsum(centres + h + 1 + GAP, centres + h + 1 + GAP + FLANK)
            E = E * ((Or + 1.0) / wn) / ((Fr + 1.0) / nf)
        sc = poisson_score(O, E, xp, ginc)
        sc = xp.round(sc, SCORE_DECIMALS)      # CPU and GPU log/gamma differ in the ~14th digit: round so that
        sc[:, xp.nonzero((nl + nr) < FLANK)[0]] = 0   # ties break identically (first maximum) on both
        i = sc.argmax(axis=1)
        s = sc[rows, i]
        up = s > best[0]
        best[0][up] = s[up]
        best[1][up] = depth[i][up]
        best[2][up] = w
        best[3][up] = O[rows, i][up]
        best[4][up] = E[rows, i][up]
    return best if xp is np else [xp.asnumpy(b) for b in best]


# ── word codes ───────────────────────────────────────────────────────────────
class Words:
    """k-mer codes and validity per promoter and start index (start j <-> depth L - j)."""
    def __init__(self, S, w, Sref):
        self.S, self.w, self.Sref = S, w, Sref
        self.S0 = S.copy()
        self._csf0 = {}                                  # unmasked sequence: flank counts
        self.code, self.valid, self.Pref = {}, {}, {}
        self.valid0, self.V0, self.V = {}, {}, {}          # unmasked validity; weighted valid starts per depth
        for k in KS:
            self.code[k], self.valid[k] = self._codes(S, k)
            self.valid0[k] = self.valid[k].copy()
            r_, j_ = np.nonzero(self.valid[k])
            self.V0[k] = np.bincount(L - j_ - 1, weights=w[r_], minlength=L)
            self.V[k] = self.V0[k].copy()
            del r_, j_
            cr, vr = self._codes(Sref, k)
            rows = np.nonzero(vr)[0]
            flat = cr[vr].astype(np.int64) * L + (L - np.nonzero(vr)[1] - 1)
            self.Pref[k] = np.bincount(flat, weights=w[rows], minlength=(4 ** k) * L
                                       ).reshape(4 ** k, L).astype(np.float32)
            del cr, vr, flat, rows

    @staticmethod
    def _codes(S, k):
        J = L - k + 1
        c = np.zeros((S.shape[0], J), dtype=np.int32)
        bad = np.zeros((S.shape[0], J), dtype=bool)
        for t in range(k):
            col = S[:, t:t + J]
            bad |= col == 4
            c = (c << 2) | (col & 3).astype(np.int32)
        return c, ~bad

    def full_scan(self, log=None, zone=(1, L)):
        """best peak (centre inside zone) of every word, all k: list of (score, k, code, depth, width)."""
        pool = []
        for k in KS:
            t0 = time.time()
            J = L - k + 1
            m = self.valid[k]
            rows = np.nonzero(m)[0]
            flat = self.code[k][m].astype(np.int64) * L + (L - np.nonzero(m)[1] - 1)
            sc = np.zeros(4 ** k)
            dd = np.zeros(4 ** k)
            ww = np.zeros(4 ** k)
            P = np.bincount(flat, weights=self.w[rows], minlength=(4 ** k) * L).reshape(4 ** k, L)
            if GPU:
                # counts and running totals are built on the CPU exactly as below (same numbers, same order of
                # addition); only the window scoring, the heavy part, runs on the GPU
                for a in range(0, 4 ** k, GPU_BATCH):      # small batches: many jobs can share one GPU
                    n = min(GPU_BATCH, 4 ** k - a)
                    cs = np.zeros((n, L + 1))
                    cs[:, 1:] = np.cumsum(P[a:a + n], axis=1)
                    csr = np.zeros_like(cs)
                    csr[:, 1:] = np.cumsum(self.Pref[k][a:a + n], axis=1)
                    try:
                        b = best_peak(cp.asarray(cs), 1, L, k, csr=cp.asarray(csr), centre_range=zone, xp=cp, ginc=cp_gammainc)
                    except cp.cuda.memory.OutOfMemoryError:  # GPU full: this batch on the CPU (identical result)
                        cp.get_default_memory_pool().free_all_blocks()
                        b = best_peak(cs, 1, L, k, csr=csr, centre_range=zone)
                    sc[a:a + n], dd[a:a + n], ww[a:a + n] = b[0], b[1], b[2]
                cp.get_default_memory_pool().free_all_blocks()
            else:
                for a in range(0, 4 ** k, 4096):
                    cs = np.zeros((min(4096, 4 ** k - a), L + 1))
                    cs[:, 1:] = np.cumsum(P[a:a + 4096], axis=1)
                    csr = np.zeros_like(cs)
                    csr[:, 1:] = np.cumsum(self.Pref[k][a:a + 4096], axis=1)
                    b = best_peak(cs, 1, L, k, csr=csr, centre_range=zone)
                    sc[a:a + cs.shape[0]], dd[a:a + cs.shape[0]], ww[a:a + cs.shape[0]] = b[0], b[1], b[2]
            top = np.argsort(-sc)[:POOL_PER_K]
            pool += [(float(sc[c]), k, int(c), int(dd[c]), int(ww[c])) for c in top if sc[c] > 0]
            if log:
                log(f"    k={k}: word scan {time.time() - t0:.0f} s")
        return pool

    def avail_cs(self, k):
        """cumulative availability (current / original weighted valid starts) by depth, shape (1, L+1)."""
        a = np.divide(self.V[k], self.V0[k], out=np.ones(L), where=self.V0[k] > 0)
        ca = np.zeros((1, L + 1))
        ca[0, 1:] = np.cumsum(a)
        return ca

    def word_peak(self, k, c, zone=(1, L)):
        hit = self.code[k] == c
        r, j = np.nonzero(hit & self.valid[k])
        prof = np.bincount(L - j - 1, weights=self.w[r], minlength=L)        # index = depth-1
        csf = self._csf0.get((k, c))                    # unmasked profile never changes: cache it
        if csf is None:
            r0, j0 = np.nonzero(hit & self.valid0[k])
            csf = np.zeros((1, L + 1))
            csf[0, 1:] = np.cumsum(np.bincount(L - j0 - 1, weights=self.w[r0], minlength=L))
            self._csf0[(k, c)] = csf
        cs = np.zeros((1, L + 1))
        cs[0, 1:] = np.cumsum(prof)
        csr = np.zeros((1, L + 1))
        csr[0, 1:] = np.cumsum(self.Pref[k][c])
        b = best_peak(cs, 1, L, k, csr=csr, centre_range=zone, csf=csf, ca=self.avail_cs(k))
        return float(b[0][0]), int(b[1][0]), int(b[2][0])

    def mask(self, inst, m):
        """mask instances: list of (promoter row, start index j, length m)."""
        for i, j in inst:
            a, b = max(0, j), min(L, j + m)
            self.S[i, a:b] = 4
            for k in KS:
                lo_, hi_ = max(0, a - k + 1), min(L - k + 1, b)
                if hi_ <= lo_:
                    continue
                was = np.nonzero(self.valid[k][i, lo_:hi_])[0] + lo_
                if len(was):
                    np.subtract.at(self.V[k], L - was - 1, self.w[i])
                self.valid[k][i, lo_:hi_] = False


# ── PWM refinement ───────────────────────────────────────────────────────────
def pwm_from(S, w, inst, m, bg):
    cnt = np.zeros((m, 4))
    for (i, j), wt in zip(inst, w):
        seg = S[i, j:j + m]
        for c, b in enumerate(seg):
            if b < 4:
                cnt[c, b] += wt
    prob = (cnt + 0.5 * 4 * bg) / (cnt.sum(axis=1, keepdims=True) + 2.0)
    return prob


def scan_pwm(S, logodds, j_lo, j_hi):
    """PWM scores for starts j in [j_lo, j_hi); masked bases -> -inf."""
    m = logodds.shape[0]
    lo_ext = np.vstack([logodds.T, np.full((1, m), -1e9)])          # (5, m)
    js = np.arange(j_lo, j_hi)
    sc = np.zeros((S.shape[0], len(js)))
    for c in range(m):                      # contiguous slice (a view), same additions in the same order
        sc += lo_ext[:, c][S[:, j_lo + c:j_hi + c]]
    return sc, js


def masked_windows(S, m, j_lo, j_hi):
    """True where the m-long window starting at j (j_lo <= j < j_hi) contains a masked/padding base (code 4)."""
    cm = np.zeros((S.shape[0], j_hi - j_lo + m), dtype=np.int32)
    np.cumsum(S[:, j_lo:j_hi + m - 1] == 4, axis=1, out=cm[:, 1:])
    return (cm[:, m:m + j_hi - j_lo] - cm[:, :j_hi - j_lo]) > 0


def refine(W, seed, bg, zone=(1, L)):
    """seed = (k, code, depth, width). Returns motif dict + instances or None."""
    S, w = W.S, W.w
    k, code, d0, w0 = seed
    m = k + 2 * EXT
    # initial instances: exact seed occurrences near the peak (one per promoter, closest)
    reach = w0 // 2 + PAD
    m_ok = (W.code[k] == code) & W.valid[k]
    r, j = np.nonzero(m_ok)
    dep = L - j
    near = np.abs(dep - d0) <= reach
    r, j = r[near], j[near]
    if len(r) == 0:
        return None
    inst = {}
    for ri, ji in zip(r, j):
        if ri not in inst or abs((L - ji) - d0) < abs((L - inst[ri]) - d0):
            inst[ri] = ji
    inst = [(ri, ji - EXT) for ri, ji in inst.items() if ji - EXT >= 0 and ji - EXT + m <= L]
    d_start = d0 + EXT                       # depth of PWM start for the seed position
    width = w0
    prev = None
    motif = None
    for it in range(8):
        prob = pwm_from(S, w[[i for i, _ in inst]], inst, m, bg)
        lo = np.log2(prob / bg[None, :])
        # region: window + flanks around the current centre (PWM start depths)
        span = max(WIDTHS) // 2 + PAD + GAP + FLANK
        dlo, dhi = max(m, d_start - span), min(L, d_start + span)
        j_lo, j_hi = L - dhi, L - dlo + 1                        # start indices
        j_hi = min(j_hi, L - m + 1)
        sc0, js = scan_pwm(W.S0, lo, j_lo, j_hi)                  # unmasked: flank counts
        sc = np.where(masked_windows(S, m, j_lo, j_hi), -1e9, sc0)  # masked cohort: same scores where untouched
        scr, _ = scan_pwm(W.Sref, lo, j_lo, j_hi)
        depths = L - js                                            # decreasing
        order = np.argsort(depths)
        sc, sc0, scr, depths = sc[:, order], sc0[:, order], scr[:, order], depths[order]   # increasing depth
        okn = ((sc > -1e8) * w[:, None]).sum(axis=0)
        ok0 = ((sc0 > -1e8) * w[:, None]).sum(axis=0)
        ca = np.zeros((1, len(okn) + 1))
        ca[0, 1:] = np.cumsum(np.divide(okn, ok0, out=np.ones(len(okn)), where=ok0 > 0))
        smax = sc.max()
        flat = sc.ravel()
        samp = flat[_RNG.integers(0, flat.size, min(flat.size, 300000))]
        samp = samp[samp > -1e8]
        cand = np.unique(np.quantile(samp, [0.990, 0.995, 0.998, 0.999, 0.9995, 0.9998])) if samp.size else []
        cand = [t for t in cand if t > 0] + [smax * f for f in (0.6, 0.7, 0.8)]
        best = None
        if len(cand):
            # all candidate thresholds scored in one call (one row each); same numbers as one call per threshold
            nt, ncl = len(cand), sc.shape[1]
            cs, csf, csr = (np.zeros((nt, ncl + 1)) for _ in range(3))
            for i, t in enumerate(cand):
                cs[i, 1:] = np.cumsum(((sc >= t) * w[:, None]).sum(axis=0))
                csf[i, 1:] = np.cumsum(((sc0 >= t) * w[:, None]).sum(axis=0))
                csr[i, 1:] = np.cumsum(((scr >= t) * w[:, None]).sum(axis=0))
            b = best_peak(cs, depths[0], depths[-1], m, csr=csr, csf=csf, ca=ca,
                          centre_range=(max(zone[0], d_start - max(WIDTHS) // 2 - PAD),
                                        min(zone[1], d_start + max(WIDTHS) // 2 + PAD)))
            for i, t in enumerate(cand):                 # first threshold reaching the maximum, as before
                if best is None or b[0][i] > best[0]:
                    best = (float(b[0][i]), int(b[1][i]), int(b[2][i]), float(b[3][i]), float(b[4][i]), t)
        if best is None or best[0] <= 0:
            return None
        score, centre, width, O, E, thr = best
        h = width // 2
        in_win = (depths >= centre - h) & (depths <= centre + h)
        sub = np.where(sc[:, in_win] >= thr, sc[:, in_win], -np.inf)
        rows = np.flatnonzero(np.isfinite(sub).any(axis=1))
        cols = sub[rows].argmax(axis=1)
        dsel = depths[in_win][cols]
        inst = [(int(ri), int(L - dd)) for ri, dd in zip(rows, dsel)]
        d_start = centre
        key = frozenset(inst)
        motif = {"prob": prob, "logodds": lo, "score": score, "d_start": centre, "width": width,
                 "O": O, "E": E, "threshold": thr, "n_inst": len(inst),
                 "w_inst": float(w[[i for i, _ in inst]].sum()) if inst else 0.0, "iterations": it + 1}
        if key == prev or not inst:
            break
        prev = key
    if not inst:
        return None
    # trim low-information end columns
    prob = pwm_from(S, w[[i for i, _ in inst]], inst, m, bg)
    ic = 2 + (prob * np.log2(prob)).sum(axis=1)
    a, b = 0, m
    while a < b and ic[a] < TRIM_BITS:
        a += 1
    while b > a and ic[b - 1] < TRIM_BITS:
        b -= 1
    if b - a < 4:
        return None
    motif["prob"] = prob[a:b]
    motif["ic"] = ic[a:b]
    motif["consensus"] = consensus(prob[a:b])
    motif["length"] = b - a
    motif["window_centre"] = motif["d_start"] - a     # centre of the scoring window
    motif["d_start"] = int(round(np.median([L - j for _, j in inst]))) - a   # median instance start
    motif["seed"] = (k, decode(code, k), d0, w0)
    inst_trim = [(i, j + a) for i, j in inst]
    return motif, inst_trim, b - a


def decode(c, k):
    return "".join(BASES[(c >> (2 * (k - 1 - t))) & 3] for t in range(k))


def consensus(prob):
    out = []
    for p in prob:
        o = np.argsort(-p)
        if p[o[0]] >= 0.6:
            out.append(BASES[o[0]])
        elif p[o[0]] + p[o[1]] >= 0.75:
            pair = "".join(sorted(BASES[o[0]] + BASES[o[1]]))
            out.append({"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M"}[pair])
        else:
            out.append("N")
    return "".join(out)


# ── discovery loop ───────────────────────────────────────────────────────────
def discover(S, w, Sref, max_motifs=MAX_MOTIFS, log=print, label="", zone=(1, L)):
    global _RNG
    _RNG = np.random.default_rng(0)        # same threshold sample in every run, whatever ran before in this process
    S = S.copy()
    W = Words(S, w, Sref)
    bgc = np.bincount(S[S < 4].ravel(), minlength=4).astype(float)
    bg = bgc / bgc.sum()
    pool = W.full_scan(log, zone)
    heap = [(-s, k, c, d, wd) for s, k, c, d, wd in pool]
    heapq.heapify(heap)
    motifs = []
    tried = 0
    while heap and len(motifs) < max_motifs:
        s, k, c, d, wd = heapq.heappop(heap)
        s2, d2, w2 = W.word_peak(k, c, zone)              # re-score on current data
        if s2 < STOP_SCORE:
            continue
        if heap and s2 < -heap[0][0] - 1e-9:
            heapq.heappush(heap, (-s2, k, c, d2, w2))
            continue
        tried += 1
        res = refine(W, (k, c, d2, w2), bg, zone)
        if res is None:
            continue
        motif, inst, m = res
        W.mask(inst, m)
        if motif["score"] < STOP_SCORE:
            continue
        motif["rank"] = len(motifs) + 1
        motifs.append(motif)
        log(f"  {label}M{motif['rank']:02d} {motif['consensus']:<18} start -{motif['d_start']:<5} "
            f"w={motif['width']:<3} O={motif['O']:.0f} E={motif['E']:.0f} x{motif['O'] / max(motif['E'], 1e-9):.2f} "
            f"score {motif['score']:.1f}  ({motif['n_inst']} instances; seed {motif['seed'][1]})")
    return motifs


def fdr_select(real, nulls):
    """FDR per window class; returns set of real motif ranks passing."""
    keep = set()
    for cls in ("sharp", "broad"):
        f = (lambda m: m["width"] <= 11) if cls == "sharp" else (lambda m: m["width"] >= 21)
        rs = sorted((m["score"] for m in real if f(m)), reverse=True)
        ns = [sorted((m["score"] for m in nl if f(m)), reverse=True) for nl in nulls]
        u = math.inf
        for t in rs:
            n_r = sum(1 for x in rs if x >= t)
            n_n = np.mean([sum(1 for x in nl if x >= t) for nl in ns]) if ns else 0
            if n_n / n_r <= FDR_MAX:
                u = t
        keep |= {m["rank"] for m in real if f(m) and m["score"] >= u}
    return keep


def _discover_job(args):
    tag, S, w, Sref, max_motifs, zone = args
    ms = discover(S, w, Sref, max_motifs, label=f"[{tag}] ", zone=zone)
    for m in ms:
        m.pop("logodds", None)
    return tag, ms


# ── bench helpers ────────────────────────────────────────────────────────────
def random_motif(rng, length=10):
    while True:
        m = "".join(rng.choice(list(BASES), length))
        gc = (m.count("G") + m.count("C")) / length
        if 0.3 <= gc <= 0.7 and max(m.count(b) for b in BASES) <= 5 and "CG" not in m:
            return m


def plant(seqs, rng, fractions, positions=(30, 90, 250, 700, 1500), jitters=(0, 15), fidelity=0.75):
    seqs = [list(s) for s in seqs]
    n = len(seqs)
    scen = []
    for p in positions:
        order = rng.permutation(n)
        used = 0
        for f in fractions:
            for j in jitters:
                m = random_motif(rng)
                cnt = int(round(f * n))
                for g in order[used:used + cnt]:
                    off = int(round(rng.normal(0, j))) if j else 0
                    off = max(-2 * j, min(2 * j, off))
                    start = L - (p + off)
                    inst = [b if rng.random() < fidelity else rng.choice([x for x in BASES if x != b]) for b in m]
                    if 0 <= start and start + len(m) <= L:
                        seqs[g][start:start + len(m)] = inst
                used += cnt
                scen.append({"motif": m, "position": p, "fraction": f, "carriers": cnt, "jitter": j})
    return ["".join(s) for s in seqs], scen


def iupac_fit(a, b):
    """number of compatible positions when IUPAC strings a, b are aligned (best shift)."""
    best = 0
    for sh in range(-len(b) + 1, len(a)):
        n = 0
        for i in range(len(a)):
            jb = i - sh
            if 0 <= jb < len(b) and set(IUPAC.get(a[i], "ACGT")) & set(IUPAC.get(b[jb], "ACGT")) \
                    and a[i] != "N" and b[jb] != "N":
                n += 1
        best = max(best, n)
    return best


def core_match(motif, cons):
    """known consensus (or its reverse complement) fits inside the motif: every
    aligned motif position is a single base, or a two-base code matching a
    two-base code in the consensus; no mismatch for cores <= 7 nt, one for longer."""
    rc = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
          "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}
    for c in (cons, "".join(rc[x] for x in reversed(cons))):
        for sh in range(0, len(motif) - len(c) + 1):
            seg = motif[sh:sh + len(c)]
            max_mm = 0 if len(c) <= 7 else 1
            mm = 0
            for a, b in zip(seg, c):
                ok = (len(IUPAC[a]) == 1 and a in IUPAC[b]) or (len(IUPAC[a]) == 2 and a == b)
                mm += not ok
            if mm <= max_mm:
                return True
    return False


_SETS = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT", ".": "ACGT"}
_RC = str.maketrans("ACGTRYKMSWBDHVN.", "TGCAYRMKSWVHDBN.")


def word_strand(motif, patterns):
    """+1 when a core word (IUPAC; '.' = any) fits inside the motif consensus on its own strand, -1 on the reverse
    strand, 0 otherwise. A motif position fits when it is a one- or two-base code sharing a base with the word."""
    def fits(c, p):
        for i in range(len(c) - len(p) + 1):
            if all(q == "." or (len(_SETS[a]) <= 2 and set(_SETS[a]) & set(_SETS[q])) for a, q in zip(c[i:i + len(p)], p)):
                return True
        return False
    for p in patterns:
        if fits(motif, p):
            return 1
    for p in patterns:
        if fits(motif.translate(_RC)[::-1], p):
            return -1
    return 0


def planted_match(cons, planted):
    """>= 7 single-base consensus positions equal to the planted motif and <= 1
    conflicting single base, at the best shift."""
    for sh in range(-len(planted) + 1, len(cons)):
        eq = bad = 0
        for i, ch in enumerate(cons):
            j = i - sh
            if 0 <= j < len(planted) and ch in "ACGT":
                eq += ch == planted[j]
                bad += ch != planted[j]
        if eq >= 7 and bad <= 1:
            return True
    return False


def similar_motif(a, b):
    """two consensuses share >= 5 single-base positions with <= 1 conflict (either strand)."""
    rc = str.maketrans("ACGTRYSWKMN", "TGCAYRSWMKN")
    for bb in (b, b.translate(rc)[::-1]):
        for sh in range(-len(bb) + 1, len(a)):
            eq = bad = 0
            for i, ch in enumerate(a):
                j = i - sh
                if 0 <= j < len(bb) and ch in "ACGT" and bb[j] in "ACGT":
                    eq += ch == bb[j]
                    bad += ch != bb[j]
            if eq >= 5 and bad <= 1:
                return True
    return False


def informative(c):
    return sum(ch in "ACGT" for ch in c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--cohort", default="C1")
    ap.add_argument("--gene-class", default=None, choices=[None, "protein_coding", "non_coding", "pseudogene"])
    ap.add_argument("--primary-only", action="store_true", help="primary anchors only (useful for C4)")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--nulls", type=int, default=5)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--max-motifs", type=int, default=MAX_MOTIFS, help="per zone")
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--zones", default="core,distal")
    ap.add_argument("--fractions", default="0.0005,0.001,0.002,0.005,0.01,0.02,0.05")
    ap.add_argument("--fidelities", default="0.75,0.85,0.95")
    ap.add_argument("--jitters", default="0,15")
    ap.add_argument("--tag", default="")
    ap.add_argument("--jobs", type=int, default=1,
                    help="discovery runs executed in parallel; ~4-5 GB RAM each")
    a = ap.parse_args()
    root = Path(a.root)
    out = root / "datasets" / "motifs"
    (out / "nulls").mkdir(parents=True, exist_ok=True)
    zones = {z: ZONES[z] for z in a.zones.split(",")}
    T0 = time.time()
    ids, seqs, w, cls = load_cohort(root, a.cohort, a.gene_class, a.primary_only)
    label = a.cohort + (f"_{a.gene_class}" if a.gene_class else "") + ("_primary" if a.primary_only else "")
    a.cohort = label                      # used in every output and cache name below
    full = sum(1 for x in seqs if not x.startswith("n"))
    print(f"{label}: {len(seqs):,} promoters ({full:,} full length), {w.sum():,.0f} independent; zones " +
          ", ".join(f"{z} -{lo}..-{hi}" for z, (lo, hi) in zones.items()))

    t0 = time.time()
    ref_seqs = null_cohort(seqs, 2999, out / "nulls" / f"null50_{a.cohort}_reference.txt.gz", a.workers)
    Sref = encode(ref_seqs)
    del ref_seqs
    print(f"Composition reference ready ({time.time() - t0:.0f} s)")

    jobs, null_motifs = [], {z: [None] * a.nulls for z in zones}
    null_S = {}
    for z, zr in zones.items():
        for r in range(a.nulls):
            cache = out / "nulls" / f"motifs_null_{a.cohort}_{z}_rep{r}_m{a.max_motifs}.json"
            if cache.exists():
                null_motifs[z][r] = json.load(open(cache))
                print(f"Null {z} {r + 1}: {len(null_motifs[z][r])} motifs (cached)")
                continue
            if r not in null_S:
                ns = null_cohort(seqs, 2000 + r, out / "nulls" / f"null50_{a.cohort}_rep{r}.txt.gz", a.workers)
                null_S[r] = encode(ns)
            jobs.append((f"{z}:null{r}", null_S[r], w, Sref, a.max_motifs, zr))
        jobs.append((f"{z}:real", encode(seqs), w, Sref, a.max_motifs, zr))
    scen_all = {}
    if a.bench:
        fr = tuple(float(x) for x in a.fractions.split(","))
        jt = tuple(int(x) for x in a.jitters.split(","))
        positions = (30, 90, 250, 700, 1500)
        for fi, fid in enumerate(float(x) for x in a.fidelities.split(",")):
            rng = np.random.default_rng(a.seed + 100 * fi)
            planted, scen = plant(seqs, rng, fr, positions=positions, jitters=jt, fidelity=fid)
            for sc in scen:
                sc["fidelity"] = fid
            scen_all[fid] = scen
            Sp = encode(planted)
            for z, zr in zones.items():
                n_in = sum(1 for sc in scen if zr[0] <= sc["position"] <= zr[1])
                jobs.append((f"{z}:planted{fid}", Sp, w, Sref, a.max_motifs + n_in, zr))
    print(f"Running {len(jobs)} discovery runs, {min(a.jobs, len(jobs))} at a time")
    results = {}
    if a.jobs > 1:
        with cf.ProcessPoolExecutor(max_workers=min(a.jobs, len(jobs))) as ex:
            for tag, ms in ex.map(_discover_job, jobs):
                results[tag] = ms
                print(f"  finished {tag}: {len(ms)} motifs ({time.time() - T0:.0f} s)")
    else:
        for j in jobs:
            tag, ms = _discover_job(j)
            results[tag] = ms
    for z in zones:
        for r in range(a.nulls):
            if null_motifs[z][r] is None:
                slim = [{k: v for k, v in m.items() if k in ("score", "width", "d_start", "consensus", "O", "E", "rank")}
                        for m in results[f"{z}:null{r}"]]
                json.dump(slim, open(out / "nulls" / f"motifs_null_{a.cohort}_{z}_rep{r}_m{a.max_motifs}.json", "w"))
                null_motifs[z][r] = slim

    def rows_of(motifs, keep, zone):
        rows = []
        for m in motifs:
            known = [nm for nm, cons, (lo, hi) in KNOWN
                     if lo - 5 <= m["d_start"] <= hi + 15 and core_match(m["consensus"], cons)]
            rows.append({"zone": zone, "rank": m["rank"], "pass_fdr": m["rank"] in keep, "consensus": m["consensus"],
                         "start_depth": m["d_start"], "length": m["length"], "window": m["width"],
                         "O": round(m["O"], 1), "E": round(m["E"], 1), "ratio": round(m["O"] / max(m["E"], 1e-9), 2),
                         "score": round(m["score"], 1), "instances": m["n_inst"],
                         "weighted_instances": round(m["w_inst"], 1), "seed": m["seed"][1], "known": ",".join(known),
                         "pwm_ACGT": ";".join(",".join(f"{x:.3f}" for x in row) for row in m["prob"])})
        return rows

    def write(rows, name):
        with open(out / name, "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["zone"], delimiter="\t")
            wr.writeheader()
            wr.writerows(rows)

    real_rows = []
    for z in zones:
        real = results[f"{z}:real"]
        real_rows += rows_of(real, fdr_select(real, null_motifs[z]), z)
    write(real_rows, f"hpsad_motifs_{a.cohort}_zones{a.tag}.tsv")
    summary = {"cohort": a.cohort, "promoters": len(seqs), "independent": round(float(w.sum()), 1),
               "zones": zones, "nulls": a.nulls,
               "null_max_scores": {z: [max([m["score"] for m in nm] or [0]) for nm in null_motifs[z]] for z in zones},
               "motifs": real_rows,
               "parameters": {"widths": WIDTHS, "flank": FLANK, "gap": GAP, "ks": KS, "ext": EXT,
                              "stop_score": STOP_SCORE, "max_motifs_per_zone": a.max_motifs, "fdr": FDR_MAX,
                              "null": "dinucleotide shuffle, 50-nt blocks, random phase"}}

    if a.bench:
        bench_rows, recall = [], []
        real_pass = [r["consensus"] for r in real_rows if r["pass_fdr"]]
        for fid, scen in scen_all.items():
            prow = []
            for z in zones:
                pm = results[f"{z}:planted{fid}"]
                prow += rows_of(pm, fdr_select(pm, null_motifs[z]), z)
            for r in prow:
                r["fidelity"] = fid
            bench_rows += prow
            for sc in scen:
                zone = next((z for z, (lo, hi) in zones.items() if lo <= sc["position"] <= hi), None)
                if zone is None:
                    continue
                hits = [r for r in prow if r["zone"] == zone and r["pass_fdr"] and planted_match(r["consensus"], sc["motif"])
                        and abs(r["start_depth"] - sc["position"]) <= 2 * sc["jitter"] + r["window"] // 2 + 8]
                sc.update(zone=zone, recovered=bool(hits), motif_ranks=[r["rank"] for r in hits])
                recall.append(sc)
            passing = [r for r in prow if r["pass_fdr"]]
            pl = [r for r in passing if any(planted_match(r["consensus"], sc["motif"]) for sc in scen)]
            rl = [r for r in passing if r not in pl and any(similar_motif(r["consensus"], c) for c in real_pass)]
            other = [r for r in passing if r not in pl and r not in rl]
            summary.setdefault("bench_false", {})[str(fid)] = {"passing": len(passing), "planted": len(pl),
                                                               "also_real": len(rl), "other": len(other),
                                                               "other_motifs": [(r["zone"], r["consensus"], r["start_depth"])
                                                                                for r in other]}
        write(bench_rows, f"hpsad_motifs_{a.cohort}_planted_zones{a.tag}.tsv")
        summary["bench"] = recall
        print("\nRECALL (recovered / planted), by zone, fidelity, carriers and jitter:")
        for z in zones:
            pos = sorted({sc["position"] for sc in recall if sc["zone"] == z})
            print(f"  {z} zone (positions {', '.join('-' + str(p) for p in pos)}):")
            fids = sorted(scen_all)
            jts = sorted({sc["jitter"] for sc in recall})
            hdr = "    carriers  " + "  ".join(f"f{f:.2f}/j{j:<2}" for f in fids for j in jts)
            print(hdr)
            for car in sorted({sc["carriers"] for sc in recall}):
                cells = []
                for f in fids:
                    for j in jts:
                        v = [sc for sc in recall if sc["zone"] == z and sc["fidelity"] == f and sc["jitter"] == j
                             and sc["carriers"] == car]
                        cells.append(f"{sum(x['recovered'] for x in v)}/{len(v)}".rjust(10))
                print(f"    {car:>8}  " + "  ".join(cells))
        print("\nPlanted runs - motifs passing FDR: planted / also in real run / other")
        for fid, d in summary["bench_false"].items():
            print(f"  fidelity {fid}: {d['passing']} = {d['planted']} / {d['also_real']} / {d['other']}  {d['other_motifs']}")

    json.dump(summary, open(out / (f"hpsad_motifs_bench_zones{a.tag}_summary.json" if a.bench
                                   else f"hpsad_motifs_{a.cohort}_zones{a.tag}_summary.json"), "w"),
              indent=1, default=float)
    print(f"\nReal data by zone (null max scores: " +
          "; ".join(f"{z} " + ",".join(f"{x:.0f}" for x in v) for z, v in summary["null_max_scores"].items()) + ")")
    for r in real_rows:
        print(f"  {r['zone']:<6} {'*' if r['pass_fdr'] else ' '} M{r['rank']:02d} {r['consensus']:<18} -{r['start_depth']:<5} "
              f"w={r['window']:<3} x{r['ratio']:<6} score {r['score']:<7} n={r['instances']:<6} {r['known']}")
    print(f"Total {time.time() - T0:.0f} s")


if __name__ == "__main__":
    main()
