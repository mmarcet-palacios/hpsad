#!/usr/bin/env python3
"""
hpsad_kmer_enrichment_v7.py  -  words enriched ANYWHERE in the promoter (no positional peak needed)

The discovery statistic compares a window with its own flanks, so an element enriched evenly along
the promoter is invisible to it. Here every 8-mer (both strands merged: a word and its reverse
complement count as one) is scored by the family-weighted number of promoters that contain it
anywhere in a window, against the same count in the ten shuffled cohorts (dinucleotide shuffles in
50-nt blocks: local composition, CpG included, is preserved). Words explained by known elements are
set aside; for the rest, the positional profile (25-nt bins) tells whether the enrichment is flat or
has a peak (peakiness = highest bin O/E / overall O/E).

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_kmer_enrichment_v7.py                 (window -300..+100, C1 with evidence)
Outputs (datasets\\profiles):
    kmer_enrichment_C1ev_down_v7.tsv    all 8-mers: O, E, sd, O/E, -log10 P, known element, peakiness
"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

import numpy as np
from scipy.special import gammainc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402

L, SHIFT, K = M.L, 100, 8
KNOWN = [("TATA box", "TATAAA"), ("TATA box", "TATAWAW"), ("CCAAT box", "CCAAT"), ("GC box", "GGGCGG"),
         ("GC box", "GGGGCGGGG"), ("CRE", "TGACG"), ("ETS", "CCGGAAG"), ("ETS", "GGAAGT"), ("ETS", "CGGAAR"),
         ("NRF1", "GCGCA"), ("NRF1", "TGCGC"), ("THAP11/Ronin", "ACTACA"), ("ZBTB33", "TCTCGCGAGA"), ("ZBTB33", "CTCGCGAG"),
         ("AP-1", "TGASTCA"), ("GGAGGAG", "GGAGGAG"), ("E-box", "CACGTG"), ("YY1", "AAGATGGC"), ("YY1", "GCCATNTT"),
         ("HNF1", "GTTAAT"), ("Initiator", "TCAKTY"), ("start codon", "CCATGG"), ("splice donor", "GTRAGT")]
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]",
         "M": "[AC]", "N": "[ACGT]"}
COMP = str.maketrans("ACGT", "TGCA")


def col(pos):
    d = SHIFT - pos + 1 if pos > 0 else SHIFT - pos
    return L - d


def word_codes(S, j0, j1):
    """8-mer codes for starts j0..j1-1 (canonical: min of word and reverse complement); -1 where masked."""
    n = j1 - j0
    code = np.zeros((S.shape[0], n), np.int64)
    rcode = np.zeros((S.shape[0], n), np.int64)
    bad = np.zeros((S.shape[0], n), bool)
    for t in range(K):
        b = S[:, j0 + t:j1 + t].astype(np.int64)
        bad |= b > 3
        b = np.minimum(b, 3)
        code = code * 4 + b
        rcode = rcode + (3 - b) * (4 ** t)
    c = np.minimum(code, rcode)
    c[bad] = -1
    return c


def presence(S, w, j0, j1):
    c = word_codes(S, j0, j1)
    rows = np.repeat(np.arange(S.shape[0]), c.shape[1])
    c = c.ravel()
    ok = c >= 0
    u = np.unique(rows[ok] * (4 ** K) + c[ok])
    return np.bincount(u % (4 ** K), weights=w[u // (4 ** K)], minlength=4 ** K)


def decode(c):
    return "".join("ACGT"[(c >> (2 * (K - 1 - t))) & 3] for t in range(K))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--window", default="-300,100")
    ap.add_argument("--nulls", type=int, default=10)
    ap.add_argument("--tag", default="_v7", help="tag of the C1ev_down run whose shuffled cohorts are used")
    a = ap.parse_args()
    root = Path(a.root)
    mdir, out = root / "datasets" / "motifs", root / "datasets" / "profiles"
    lo, hi = (int(x) for x in a.window.split(","))
    ids, seqs, w, _ = M.load_cohort(root, "C1", "!pseudogene", False, sources=["EPDnew", "CAGE"],
                                    segments="hpsad_segments_v7_down100.tsv.gz")
    S = M.encode(seqs)
    shuf = [M.encode(gzip.open(mdir / "nulls" / f"null50_C1ev_down{a.tag}_rep{r}.txt.gz", "rt").read().split("\n"))
            for r in range(a.nulls)]
    j0, j1 = col(lo), col(hi) - K + 2                  # word starts so that the whole word lies in the window
    O = presence(S, w, j0, j1)
    En = np.array([presence(Z, w, j0, j1) for Z in shuf])
    E, sd = En.mean(0), En.std(0)
    canon = np.array([c for c in range(4 ** K) if c == min(c, int("".join(str(3 - int(x)) for x in np.base_repr(c, 4).zfill(K))[::-1], 4))])
    words = [decode(c) for c in canon]
    known = []
    for wd in words:
        rcw = wd.translate(COMP)[::-1]
        k = ""
        for nm, cons in KNOWN:
            rx = "".join(IUPAC[ch] for ch in cons)
            rrc = "".join(IUPAC[ch] for ch in cons.translate(str.maketrans("ACGTRYKM", "TGCAYRMK"))[::-1])
            if re.search(rx, wd) or re.search(rx, rcw) or re.search(rrc, wd):
                k = nm
                break
        known.append(k)
    o, e = O[canon], E[canon]
    sc = np.where(o > e, -np.log10(np.maximum(gammainc(np.maximum(o, 1e-9), np.maximum(e, 1e-9)), 1e-300)), 0.0)
    ratio = (o + 0.5) / (e + 0.5)
    # peakiness for candidate words: positional O/E in 25-nt bins
    sig = np.nonzero((ratio >= 1.5) & (sc >= 10))[0]
    sig = sig[np.argsort(-sc[sig])]
    cand = [i for i in sig if not known[i]][:60] + [i for i in sig if known[i]][:20]     # peakiness for the top words
    bins = list(range(lo, hi, 25))
    cw = word_codes(S, j0, j1)
    cwn = [word_codes(Z, j0, j1) for Z in shuf]
    starts = np.arange(j0, j1)
    pos = np.array([-(L - j - SHIFT) if L - j > SHIFT else SHIFT + 1 - (L - j) for j in starts])
    peak = np.full(len(canon), np.nan)
    best_bin = np.full(len(canon), 0)
    for i in cand:
        c = canon[i]
        hr = (cw == c)
        hn = [(x == c) for x in cwn]
        vals = []
        for b0 in bins:
            m = (pos >= b0) & (pos < b0 + 25)
            ob = w @ hr[:, m].any(1)
            eb = np.mean([w @ h[:, m].any(1) for h in hn])
            vals.append((ob + 0.5) / (eb + 0.5))
        peak[i] = max(vals) / ratio[i]
        best_bin[i] = bins[int(np.argmax(vals))]
    rows = []
    for i in np.argsort(-sc):
        rows.append(dict(word=words[i], O=round(float(o[i]), 1), E=round(float(e[i]), 1), E_sd=round(float(sd[canon[i]]), 1),
                         ratio=round(float(ratio[i]), 2), score=round(float(sc[i]), 1), known=known[i],
                         peakiness=("" if np.isnan(peak[i]) else round(float(peak[i]), 2)),
                         best_bin=("" if np.isnan(peak[i]) else f"{best_bin[i]:+d}..{best_bin[i] + 24:+d}")))
    out.mkdir(parents=True, exist_ok=True)
    with open(out / f"kmer_enrichment_C1ev_down{a.tag}.tsv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    ok = [r for r in rows if r["ratio"] >= 1.5 and r["score"] >= 10]
    print(f"{len(ids):,} promoters, window {lo:+d}..{hi:+d}; {len(canon):,} words; "
          f"{len(ok)} enriched >= 1.5-fold with -log10 P >= 10, {sum(1 for r in ok if not r['known'])} not explained by known elements")
    print(f"{'word':<10}{'O':>8}{'E':>8}{'O/E':>7}{'score':>8}{'peaky':>7}  {'best bin':<12}known")
    for r in [r for r in ok if not r["known"]][:40]:
        print(f"{r['word']:<10}{r['O']:>8}{r['E']:>8}{r['ratio']:>7}{r['score']:>8}{r['peakiness']:>7}  {r['best_bin']:<12}")
    print("\nknown-element words, top 15:")
    for r in [r for r in ok if r["known"]][:15]:
        print(f"{r['word']:<10}{r['O']:>8}{r['E']:>8}{r['ratio']:>7}{r['score']:>8}{r['peakiness']:>7}  {r['best_bin']:<12}{r['known']}")


if __name__ == "__main__":
    main()
