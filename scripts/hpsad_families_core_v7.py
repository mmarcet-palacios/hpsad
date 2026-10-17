#!/usr/bin/env python3
"""
hpsad_families_core_v7.py  -  families extended to copies shared only in the discovery window

hpsad_families_v7.py links segments that are copies over >= 500 nt of the upstream 2,000 nt. Once the
sequence is extended to +100 (extract_downstream_v7.py), multicopy genes whose UPSTREAM sequences
differ but whose first 100 nt (the gene body: tRNA, snRNA, snoRNA ... copies) are identical stay
unlinked, and each copy counts as an independent promoter. This script adds those links.

Window: the last 400 nt of each composite sequence (-300..-1, +1..+100), upper case only (repeat-
masked bases never count). Words: 16-mers with at least 3 different bases and no base above 75%
(low-complexity runs are not evidence of a copy). Unlike the upstream rule, words shared by many
segments are NOT discarded: a 400-copy tRNA family is exactly what has to be found.
Pairs whose windows overlap in the genome (alternative anchors of one gene) are not copies and are not linked.
Link: two segments share >= MIN_WORDS 16-mers on one alignment diagonal (offset tolerance +-DIAG_TOL),
spread over >= MIN_SPAN nt. Families = connected groups of the union of these links and the upstream
families of hpsad_families_v7.tsv; weight = 1 / family size.

Inputs : datasets/cohorts/hpsad_segments_v7_down100.tsv.gz, datasets/families/hpsad_families_v7.tsv
Outputs: datasets/families/hpsad_families_v7_down100.tsv          one row per segment (same columns + core_links)
         datasets/families/hpsad_families_v7_down100_summary.json  new links, largest new families, by biotype
Usage  : python scripts\\hpsad_families_core_v7.py
"""
import argparse
import csv
import gzip
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

K, WIN = 16, 400
MIN_WORDS, MIN_SPAN, DIAG_TOL = 5, 50, 6
MAX_DF = 5000                     # guard against pathological words only


def words(s):
    """(code, pos) of acceptable 16-mers in s (upper-case ACGT only)."""
    out = []
    for p in range(len(s) - K + 1):
        w = s[p:p + K]
        if not w.isupper() or any(c not in "ACGT" for c in w):
            continue
        cnt = Counter(w)
        if len(cnt) < 3 or max(cnt.values()) > 0.75 * K:
            continue
        out.append((w, p))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    root = Path(a.root)
    t0 = time.time()
    csv.field_size_limit(10 ** 8)
    ids, bio, win, loc, span = [], [], [], [], []
    with gzip.open(root / "datasets" / "cohorts" / "hpsad_segments_v7_down100.tsv.gz", "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ids.append(r["anchor_id"])
            bio.append(r["biotype"])
            win.append(r["sequence_softmasked"][-WIN:])
            t = int(r["anchor_tss"])
            lo_, hi_ = (t - (WIN - 100), t + 99) if r["strand"] == "+" else (t - 99, t + (WIN - 100))
            loc.append((r["chrom"], lo_, hi_))
            span.append((r["chrom"], int(r["seg_start"]) - 100, int(r["seg_end"]) + 100))
    n = len(ids)
    idx = {x: i for i, x in enumerate(ids)}
    print(f"{n:,} segments; window = last {WIN} nt ({time.time() - t0:.0f} s)")

    occ = defaultdict(list)                       # word -> [(seg, pos)]
    for i, s in enumerate(win):
        seen = set()
        for w, p in words(s):
            if w not in seen:                     # first occurrence per segment
                seen.add(w)
                occ[w].append((i, p))
    print(f"{len(occ):,} distinct words ({time.time() - t0:.0f} s)")

    # candidate pairs: count shared words per (i, j, diagonal bin)
    keys = []
    for w, lst in occ.items():
        m = len(lst)
        if m < 2 or m > MAX_DF:
            continue
        arr = np.array(lst, dtype=np.int64)
        ii, jj = np.triu_indices(m, 1)
        a_, b_ = arr[ii], arr[jj]
        lo = np.minimum(a_[:, 0], b_[:, 0])
        hi = np.maximum(a_[:, 0], b_[:, 0])
        sw = a_[:, 0] > b_[:, 0]
        pa = np.where(sw, b_[:, 1], a_[:, 1])
        pb = np.where(sw, a_[:, 1], b_[:, 1])
        d = (pb - pa + WIN) // (2 * DIAG_TOL)
        keys.append(np.stack([lo, hi, d, pa], axis=1))
    K_ = np.concatenate(keys)
    print(f"{len(K_):,} shared-word hits ({time.time() - t0:.0f} s)")
    order = np.lexsort((K_[:, 3], K_[:, 2], K_[:, 1], K_[:, 0]))
    K_ = K_[order]
    grp = np.r_[True, np.any(K_[1:, :3] != K_[:-1, :3], axis=1)]
    gid = np.cumsum(grp) - 1
    cnt = np.bincount(gid)
    pmin = np.minimum.reduceat(K_[:, 3], np.nonzero(grp)[0])
    pmax = np.maximum.reduceat(K_[:, 3], np.nonzero(grp)[0])
    first = K_[grp]
    ok = (cnt >= MIN_WORDS) & (pmax - pmin + K >= MIN_SPAN)
    links = {(int(x), int(y)) for x, y in first[ok][:, :2]}
    same = {(x, y) for x, y in links if loc[x][0] == loc[y][0] and loc[x][1] <= loc[y][2] and loc[y][1] <= loc[x][2]}
    links -= same                                  # overlapping windows of the same locus are the same DNA, not copies
    print(f"{len(same):,} pairs dropped as the same locus (overlapping windows)")
    print(f"{len(links):,} linked pairs ({time.time() - t0:.0f} s)")

    # union with upstream families
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    fam_rows = {}
    with open(root / "datasets" / "families" / "hpsad_families_v7.tsv", encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        byfam = defaultdict(list)
        for r in rd:
            fam_rows[r["anchor_id"]] = r
            if r["anchor_id"] in idx:
                byfam[r["family_id"]].append(idx[r["anchor_id"]])
    for mem in byfam.values():
        for x in mem[1:]:
            parent[find(x)] = find(mem[0])
    old_root = [find(i) for i in range(n)]
    new_links = 0
    for x, y in links:
        if find(x) != find(y):
            new_links += 1
            parent[find(x)] = find(y)
    comp = defaultdict(list)
    for i in range(n):
        comp[find(i)].append(i)
    def n_loci(members):                            # overlapping windows of one locus count once
        iv = sorted(span[i] for i in members)
        k, last = 0, None
        for c, lo, hi in iv:
            if last and c == last[0] and lo <= last[2]:
                last = (c, last[1], max(hi, last[2]))
            else:
                k += 1
                last = (c, lo, hi)
        return k
    size = {r: n_loci(m) for r, m in comp.items()}
    core_deg = Counter()
    for x, y in links:
        core_deg[x] += 1
        core_deg[y] += 1

    out = root / "datasets" / "families"
    names = {r: f"G{k:05d}" for k, r in enumerate(sorted(comp, key=lambda r: -size[r]), 1)}
    changed = 0
    with open(out / "hpsad_families_v7_down100.tsv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols + ["family_down100", "family_size_down100", "core_links"], delimiter="\t")
        wr.writeheader()
        for i, aid in enumerate(ids):
            r = dict(fam_rows.get(aid, {"anchor_id": aid}))
            rt = find(i)
            r["family_down100"] = names[rt]
            r["family_size_down100"] = size[rt]
            r["core_links"] = core_deg[i]
            old_w = float(r.get("weight") or 1)
            r["weight"] = round(1 / size[rt], 6)
            changed += abs(old_w - r["weight"]) > 1e-9
            wr.writerow(r)
    grown = sorted((m for m in comp.values() if len(m) > 1 and len({old_root[i] for i in m}) > 1), key=len, reverse=True)
    summary = {"segments": n, "core_linked_pairs": len(links), "links_joining_families": new_links,
               "segments_with_changed_weight": changed,
               "independent_before": round(sum(float(fam_rows[x]["weight"]) for x in ids if x in fam_rows), 1),
               "independent_after": round(sum(1 / size[find(i)] for i in range(n)), 1),
               "segments_core_linked_by_biotype": dict(Counter(bio[i] for i in core_deg).most_common()),
               "largest_grown_families": [{"size": len(m), "biotypes": dict(Counter(bio[i] for i in m).most_common(4)),
                                           "examples": [ids[i].split("|")[0].replace("gene-", "") for i in m[:6]]}
                                          for m in grown[:15]]}
    json.dump(summary, open(out / "hpsad_families_v7_down100_summary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "largest_grown_families"}, indent=1))
    for g in summary["largest_grown_families"]:
        print(f"  {g['size']:>4}  {g['biotypes']}  {', '.join(g['examples'])}")
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
