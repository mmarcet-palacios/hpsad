#!/usr/bin/env python3
"""
hpsad_gpu_check.py  -  does the GPU word scan give the same answer as the CPU, and how much faster?

Runs the full word scan (k = 6, 7, 8) for both zones on a subset of C1_evidence promoters, once on
the GPU and once on the CPU, and compares the candidate pools (same words, same best depth and
width, scores equal to within 1e-6 relative) and the first 5 motifs discovered per zone. Prints the timings.

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_gpu_check.py                 (3,000 promoters)
    python scripts\\hpsad_gpu_check.py --n 12000       (the whole cohort)
"""
import argparse
import os
import sys
import time
from pathlib import Path

os.environ["HPSAD_GPU"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--n", type=int, default=3000)
    a = ap.parse_args()
    if not M.GPU:
        sys.exit("GPU backend not available (see the message above); install CuPy first.")
    ids, seqs, w, cls = M.load_cohort(Path(a.root), "C1", "!pseudogene", False, sources=["EPDnew", "CAGE"])
    seqs, w = seqs[:a.n], w[:a.n]
    S = M.encode(seqs)
    Sref = M.encode(M._shuffle_chunk((seqs, 7)))
    W = M.Words(S, w, Sref)
    M.GPU = True                                   # warm-up: first GPU call compiles its kernels
    cs = M.np.zeros((8, M.L + 1))
    M.best_peak(M.cp.asarray(cs), 1, M.L, 6, csr=M.cp.asarray(cs), xp=M.cp, ginc=M.cp_gammainc)
    ok = True
    for zname, zone in M.ZONES.items():
        res = {}
        for mode in ("gpu", "cpu"):
            M.GPU = mode == "gpu"
            t = time.time()
            res[mode] = (W.full_scan(zone=zone), time.time() - t)
        M.GPU = True
        g, c = res["gpu"][0], res["cpu"][0]
        gd = {(k, code): (s, d, wd) for s, k, code, d, wd in g}
        cd = {(k, code): (s, d, wd) for s, k, code, d, wd in c}
        common = set(gd) & set(cd)
        same_pos = sum(gd[x][1:] == cd[x][1:] for x in common)
        rel = max((abs(gd[x][0] - cd[x][0]) / max(abs(cd[x][0]), 1e-9) for x in common), default=0.0)
        top_g = [x for x in sorted(gd, key=lambda x: -gd[x][0])[:200]]
        top_c = [x for x in sorted(cd, key=lambda x: -cd[x][0])[:200]]
        print(f"{zname:<7} words: GPU {len(g)}, CPU {len(c)}, shared {len(common)}; same depth+width {same_pos}/{len(common)}; "
              f"max relative score difference {rel:.2e}; top-200 identical: {top_g == top_c}")
        print(f"        time: GPU {res['gpu'][1]:.1f} s, CPU {res['cpu'][1]:.1f} s  ->  {res['cpu'][1] / max(res['gpu'][1], 1e-9):.1f}x")
        ok &= len(common) == len(cd) == len(gd) and same_pos == len(common) and rel < 1e-6 and top_g == top_c
    # what matters in the end: the motifs discovered (5 per zone)
    for zname, zone in M.ZONES.items():
        found = {}
        for mode in ("gpu", "cpu"):
            M.GPU = mode == "gpu"
            t = time.time()
            ms = M.discover(S, w, Sref, 5, log=lambda *x, **y: None, label="check", zone=zone)
            found[mode] = ([(m["consensus"], m["d_start"], round(m["score"], 4)) for m in ms], time.time() - t)
        M.GPU = True
        same = found["gpu"][0] == found["cpu"][0]
        print(f"{zname:<7} discovery (5 motifs): identical {same}; GPU {found['gpu'][1]:.0f} s, CPU {found['cpu'][1]:.0f} s")
        if not same:
            print("   GPU:", found["gpu"][0])
            print("   CPU:", found["cpu"][0])
        ok &= same
    print("\nRESULT:", "GPU matches CPU - safe to use --gpu" if ok else "GPU DIFFERS from CPU - do not use --gpu yet; paste this output")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
