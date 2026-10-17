#!/usr/bin/env python3
"""
hpsad_families_v7.py  -  HPSAD step 3: promoter families (shared ancestry)

QUESTION: are two promoter segments copies of each other (duplicated genes,
pseudogenes copied from a parent gene, tandem arrays)? If so they are one
family and together count as ONE observation in discovery.

METHOD (plagiarism-check logic)
  1. Words: every 16-nt word (16-mer) in each segment. A 16-mer is long enough
     that two unrelated segments almost never share one by chance
     (4^16 = 4.3 billion possible words).
  2. Ignore repeat words: 16-mers overlapping any repeat-masked base (lower case
     in the NCBI FASTA) are never used, and any remaining 16-mer present in more
     than MAX_DF segments is also ignored. Interspersed repeats (Alu, L1, ERV ...)
     are shared by unrelated promoters, so they are not evidence of a copy -
     like ignoring "the" and "and" when checking two essays for plagiarism.
  3. Stage 1 (fast, permissive): pairs sharing >= MIN_SHARED rare words and
        s = shared rare 16-mers / min(rare 16-mers in A, rare 16-mers in B) >= S_CAND
  4. Stage 2 (verification): a true copy shares its words IN THE SAME ORDER
     along the segment (one alignment diagonal), spread over >= MIN_SPAN nt.
     A shared repeat or chance match gives a short patch or scattered words.
     Identity over the spanned region is estimated from the density of exact
     shared 16-mers (density ~ p^16); pairs with identity >= MIN_IDENT are linked.
  5. Families = connected groups of linked segments. weight = 1 / family size.

Thresholds are calibrated by --selftest: real promoters are copied with known
identity (random substitutions plus small insertions/deletions) and the
detection rate is measured, together with false links to shuffled copies.

Inputs : datasets/cohorts/hpsad_segments_v7.tsv.gz
Outputs: datasets/families/hpsad_families_v7.tsv   one row per segment
         datasets/families/hpsad_family_pairs_v7.tsv linked pairs with scores
         datasets/families/hpsad_families_v7_summary.json
Usage  : python scripts\\hpsad_families_v7.py              (full run)
         python scripts\\hpsad_families_v7.py --selftest   (calibration only)
"""

import argparse
import csv
import gzip
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

K = 16
MAX_DF = 50          # words in more segments than this are treated as repeat words
S_CAND = 0.02        # stage 1: candidate pairs (very permissive)
MIN_SHARED = 20
DIAG_TOL = 60        # stage 2: shared words must lie on one diagonal (+/- indel drift)
MIN_SPAN = 500       # ... spread over >= 500 nt of both segments (or 50% of a shorter one)
MIN_COLINEAR = 0.8   # ... and >= 80% of shared words must be on that diagonal
MIN_IDENT = 0.80     # ... with estimated identity >= 80% over the spanned region
LUT = np.full(256, 255, dtype=np.uint8)
for i, ch in enumerate(b"ACGT"):
    LUT[ch] = i


def kmers(seq, with_pos=False):
    """unique 16-mer codes (uint32) of an ACGT string (and first positions)."""
    a = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(a) - K + 1
    if n <= 0:
        e = np.empty(0, dtype=np.uint32)
        return (e, np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)) if with_pos else e
    a = a.astype(np.uint32)
    code = np.zeros(n, dtype=np.uint32)
    for j in range(K):
        code = (code << 2) | a[j:j + n]
    bad = np.convolve((a == 255).astype(np.int32), np.ones(K, dtype=np.int32), "valid") > 0
    if with_pos:
        pos = np.flatnonzero(~bad)
        u, first = np.unique(code[pos], return_index=True)
        return u, pos[first], pos
    return np.unique(code[~bad])


def link_pairs(seqs, max_df=MAX_DF, s_min=S_CAND, min_shared=MIN_SHARED, log=print, verify=True):
    t0 = time.time()
    per = [kmers(s) for s in seqs]
    lens = np.array([len(x) for x in per])
    codes = np.concatenate(per)
    owner = np.repeat(np.arange(len(per), dtype=np.int32), lens)
    del per
    order = np.argsort(codes, kind="stable")
    codes, owner = codes[order], owner[order]
    del order
    starts = np.flatnonzero(np.r_[True, codes[1:] != codes[:-1]])
    df = np.diff(np.r_[starts, len(codes)])
    log(f"  {len(codes):,} word occurrences, {len(starts):,} distinct words ({time.time() - t0:.0f} s)")

    # rare-word count per segment (denominator)
    rare_word = np.repeat(df <= max_df, df)
    rare_n = np.bincount(owner[rare_word], minlength=len(seqs))

    # pairs from words shared by 2..max_df segments
    keep = (df >= 2) & (df <= max_df)
    pair_counts = Counter()
    for d in range(2, max_df + 1):
        sel = starts[keep & (df == d)]
        if len(sel) == 0:
            continue
        mem = owner[sel[:, None] + np.arange(d)]            # (words, d)
        ii, jj = np.triu_indices(d, 1)
        a, b = mem[:, ii].ravel(), mem[:, jj].ravel()
        lo, hi = np.minimum(a, b).astype(np.int64), np.maximum(a, b).astype(np.int64)
        key = lo * len(seqs) + hi
        u, c = np.unique(key, return_counts=True)
        pair_counts.update(dict(zip(u.tolist(), c.tolist())))
    log(f"  {len(pair_counts):,} candidate pairs ({time.time() - t0:.0f} s)")

    cands = []
    n = len(seqs)
    for key, shared in pair_counts.items():
        if shared < min_shared:
            continue
        i, j = divmod(key, n)
        denom = min(rare_n[i], rare_n[j])
        s = shared / denom if denom else 0.0
        if s >= s_min:
            cands.append((i, j, shared, s))
    log(f"  {len(cands):,} candidate pairs at s >= {s_min} ({time.time() - t0:.0f} s)")
    if not verify:
        return [(i, j, sh, round(s, 4), None, None) for i, j, sh, s in cands], rare_n
    links = verify_pairs(seqs, cands)
    log(f"  {len(links):,} verified links ({time.time() - t0:.0f} s)")
    return links, rare_n


def verify_pairs(seqs, cands):
    """Stage 2: shared 16-mers must be colinear (one diagonal) and spread out.
    Returns (i, j, shared, s, est_identity, span)."""
    cache = {}

    def kp(i):
        if i not in cache:
            if len(cache) > 20000:
                cache.clear()
            cache[i] = kmers(seqs[i], with_pos=True)
        return cache[i]
    out = []
    for i, j, shared, s in cands:
        ci, pi, vi = kp(i)
        cj, pj, vj = kp(j)
        _, xi, xj = np.intersect1d(ci, cj, assume_unique=True, return_indices=True)
        a, b = pi[xi], pj[xj]
        d = a - b
        med = np.median(d)
        on = np.abs(d - med) <= DIAG_TOL
        if on.sum() < MIN_SHARED or on.mean() < MIN_COLINEAR:
            continue
        a_on, b_on = a[on], b[on]
        span_a = a_on.max() - a_on.min() + K
        span_b = b_on.max() - b_on.min() + K
        span = min(span_a, span_b)
        shorter = min(len(seqs[i]), len(seqs[j]))
        if span < min(MIN_SPAN, 0.5 * shorter):
            continue
        # identity: fraction of start positions in the spanned region that begin
        # an exact shared 16-mer ~ p^16
        valid_i = np.count_nonzero((vi >= a_on.min()) & (vi <= a_on.max()))
        valid_j = np.count_nonzero((vj >= b_on.min()) & (vj <= b_on.max()))
        frac = on.sum() / max(min(valid_i, valid_j), 1)
        ident = min(frac, 1.0) ** (1 / K)
        if ident >= MIN_IDENT:
            out.append((i, j, int(shared), round(s, 4), round(ident, 4), int(span)))
    return out


def locus_counts(rows, root_of, pad=100):
    """number of distinct genomic loci per family: members whose windows overlap (e.g. the main and an alternative
    start site of one gene) are the same DNA, not copies, and count as one locus."""
    members = defaultdict(list)
    for i, r in enumerate(rows):
        members[root_of[i]].append((r["chrom"], int(r["seg_start"]) - pad, int(r["seg_end"]) + pad))
    nloc = {}
    for root, iv in members.items():
        iv.sort()
        k, last = 0, None
        for c, lo, hi in iv:
            if last and c == last[0] and lo <= last[2]:
                last = (c, last[1], max(hi, last[2]))
            else:
                k += 1
                last = (c, lo, hi)
        nloc[root] = k
    return nloc


def families(n, links):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, j, *_ in links:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    root = [find(i) for i in range(n)]
    size = Counter(root)
    return root, size


# ── calibration ──────────────────────────────────────────────────────────────
def mutate(seq, identity, rng):
    """copy with (1-identity) of positions changed: 85% substitutions, 15% 1-3 nt indels."""
    s = list(seq)
    n_events = int(round(len(s) * (1 - identity)))
    for _ in range(n_events):
        p = int(rng.integers(0, len(s)))
        r = rng.random()
        if r < 0.85:
            new = "ACGT"[(("ACGT".index(s[p].upper()) + int(rng.integers(1, 4))) % 4)]
            s[p] = new if s[p].isupper() else new.lower()
        elif r < 0.925:
            del s[p:p + int(rng.integers(1, 4))]
        else:
            s[p:p] = list(rng.choice(list("ACGT"), int(rng.integers(1, 4))))
    return "".join(s)


def selftest(seqs, rng, log=print):
    base_idx = rng.choice(len(seqs), 300, replace=False)
    levels = [1.0, 0.95, 0.90, 0.85, 0.82, 0.80, 0.75, 0.70]
    test, origin = list(seqs), {}
    for b in base_idx:
        for lv in levels:
            origin[len(test)] = (int(b), lv)
            test.append(mutate(seqs[b], lv, rng))
        origin[len(test)] = (int(b), "shuffled")
        test.append("".join(rng.permutation(list(seqs[b]))))
        # unrelated promoter carrying a 300-nt exact chunk of b (mimics a shared repeat)
        other = seqs[int(rng.integers(0, len(seqs)))]
        st = int(rng.integers(0, max(1, len(seqs[b]) - 300)))
        origin[len(test)] = (int(b), "shared_300nt_chunk")
        test.append(other[:800] + seqs[b][st:st + 300] + other[1100:])
    log(f"Self-test: {len(base_idx)} promoters x {len(levels)} identities + shuffled copies")
    links, _ = link_pairs(test, log=log)
    found = defaultdict(set)
    for i, j, *_ in links:
        for x, y in ((i, j), (j, i)):
            if x in origin and origin[x][0] == y:
                found[origin[x][1]].add(x)
    res = {str(lv): round(len(found[lv]) / len(base_idx), 3) for lv in levels + ["shuffled", "shared_300nt_chunk"]}
    log("Detection rate of planted copies by true identity: " + json.dumps(res))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    root = Path(a.root)
    out = root / "datasets" / "families"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    csv.field_size_limit(10 ** 7)
    rows = []
    with gzip.open(root / "datasets" / "cohorts" / "hpsad_segments_v7.tsv.gz", "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["has_N"] == "0":
                rows.append(r)
    col = "sequence_softmasked" if "sequence_softmasked" in rows[0] else "sequence_5to3"
    print(f"Using column {col}" + ("" if col == "sequence_softmasked" else
          "  (WARNING: no repeat mask - re-run extract_cohort_sequences_v7.py)"))
    seqs = [r[col] for r in rows]
    for r in rows:
        r.pop("sequence_5to3", None)
        r.pop("sequence_softmasked", None)
    print(f"{len(seqs):,} N-free segments")

    c1 = [s for r, s in zip(rows, seqs) if r["cohort"] == "C1"]
    calib = selftest(c1, rng)
    if a.selftest:
        return

    print("Full run")
    links, rare_n = link_pairs(seqs)
    root_of, size_members = families(len(seqs), links)
    size = locus_counts(rows, root_of)            # family size = distinct loci (weight = 1 / loci)
    print(f"  families with >= 2 members: {sum(1 for v in size_members.values() if v > 1):,}; "
          f"with >= 2 distinct loci: {sum(1 for v in size.values() if v > 1):,}")
    fam_id = {r: f"F{k:05d}" for k, r in enumerate(sorted(set(root_of), key=lambda x: (-size[x], x)))}
    best = {}
    for i, j, shared, s, ident, span in links:
        for x, y in ((i, j), (j, i)):
            if x not in best or ident > best[x][2]:
                best[x] = (y, s, ident)

    with open(out / "hpsad_families_v7.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["anchor_id", "gene_id", "symbol", "cohort", "gene_class", "family_id",
                    "family_size", "weight", "best_partner", "best_s", "est_identity", "rare_16mers"])
        for i, r in enumerate(rows):
            sz = size[root_of[i]]
            b = best.get(i)
            w.writerow([r["anchor_id"], r["gene_id"], r["symbol"], r["cohort"], r["gene_class"],
                        fam_id[root_of[i]], sz, round(1 / sz, 6),
                        rows[b[0]]["anchor_id"] if b else "", b[1] if b else "", b[2] if b else "",
                        int(rare_n[i])])
    with open(out / "hpsad_family_pairs_v7.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["anchor_a", "anchor_b", "shared_16mers", "s", "est_identity", "aligned_span"])
        for i, j, shared, s, ident, span in links:
            w.writerow([rows[i]["anchor_id"], rows[j]["anchor_id"], shared, s, ident, span])

    # summary (family size = distinct loci)
    fam_members = defaultdict(list)
    for i, r in enumerate(rows):
        fam_members[root_of[i]].append(r)
    multi = [(rt, m) for rt, m in fam_members.items() if size[rt] > 1]
    eff = defaultdict(lambda: [0, 0.0])
    for i, r in enumerate(rows):
        k = f'{r["cohort"]}|{r["gene_class"]}'
        eff[k][0] += 1
        eff[k][1] += 1 / size[root_of[i]]
    top = sorted(multi, key=lambda x: -size[x[0]])[:25]
    summary = {
        "parameters": {"K": K, "MAX_DF": MAX_DF, "S_CAND": S_CAND, "MIN_SHARED": MIN_SHARED, "DIAG_TOL": DIAG_TOL, "MIN_SPAN": MIN_SPAN, "MIN_COLINEAR": MIN_COLINEAR, "MIN_IDENT": MIN_IDENT, "seed": a.seed,
                       "family_size": "distinct loci (overlapping windows of one locus count once)"},
        "calibration_detection_by_identity": calib,
        "segments": len(rows),
        "links": len(links),
        "families_total": len(fam_members),
        "linked_groups_with_2plus_members": sum(1 for m in fam_members.values() if len(m) > 1),
        "segments_in_multi_member_families": sum(len(m) for _, m in multi),
        "multi_member_families": len(multi),
        "family_size_distribution": dict(sorted(Counter(size[rt] for rt, _ in multi).items())),
        "segments_vs_effective_N_by_cohort_class": {k: [v[0], round(v[1], 1)] for k, v in sorted(eff.items())},
        "largest_families": [{"size": size[rt], "segments": len(m), "example_symbols": sorted({x["symbol"] for x in m})[:12],
                              "cohorts": dict(Counter(x["cohort"] for x in m))} for rt, m in top],
    }
    with open(out / "hpsad_families_v7_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps({k: summary[k] for k in ("links", "families_total", "multi_member_families",
                                              "segments_in_multi_member_families",
                                              "family_size_distribution")}, indent=1))


if __name__ == "__main__":
    main()
