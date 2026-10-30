#!/usr/bin/env python3
"""
hpsad_s1_table_v7.py  -  S1 Table: every main-cohort promoter, its sequence and every element call, with genome links

Re-runs the element calls of hpsad_architecture_v7.py (same elements, motifs, thresholds and positional zones) on the
main cohort (C1, EPDnew or CAGE start site, no pseudogenes) and records each individual match: element, motif, position
relative to the start site, strand, matched sequence, genomic coordinates (GRCh38) and a UCSC Genome Browser link.
The promoter-level presence/absence calls are checked against architecture_calls_C1ev_v7.tsv before anything is written.

Positions: +1 is the start-site base, -1 the base immediately upstream; there is no position 0.
Genomic coordinates are 1-based, inclusive, on GRCh38 (hg38).

Output: datasets/profiles/S1_Table_HPSAD_promoters.xlsx (and S1_Table_elements.tsv, S1_Table_promoters.tsv)
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
import hpsad_architecture_v7 as A                             # noqa: E402
from hpsad_syntax_v7 import cpg_class                         # noqa: E402

L, SHIFT = M.L, A.SHIFT
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
SUPPORT = {"TATA": "6 of 6 replicates", "CCAAT": "6 of 6 replicates", "GC box": "6 of 6 replicates", "ETS": "6 of 6 replicates",
           "CRE": "6 of 6 replicates", "THAP11": "6 of 6 replicates", "AP-1": "4 of 6 replicates",
           "Inr": "start-site run", "TCT": "start-site run", "YY1": "start-site run"}   # paper, Fig 5B and Fig 6


def revcomp(s):
    return s.translate(COMP)[::-1]


def chrom_name(acc):
    """RefSeq GRCh38 accession -> UCSC chromosome name (NC_000001.11 -> chr1, ..23 -> chrX, ..24 -> chrY)."""
    if acc.startswith("NC_0000"):
        n = int(acc.split(".")[0][-2:])
        return "chr" + {23: "X", 24: "Y"}.get(n, str(n))
    if acc.startswith("NC_012920"):
        return "chrM"
    return acc


def offset(c):
    """composite column (0..L-1, last column = +100) -> 0-based transcript offset from the +1 base."""
    return c - (L - SHIFT)


def rel_pos(t):
    """0-based transcript offset -> position label (+1 = start site, -1 = upstream base, no 0)."""
    return t + 1 if t >= 0 else t


def genomic(tss, strand, t):
    return tss + t if strand == "+" else tss - t


def ucsc(chrom, a, b, pad=0):
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}:{max(1, a - pad)}-{b + pad}"


def pwm_detail(S, P, thr, lo, hi, both):
    """all matches starting at positions lo..hi: list per promoter of (column, orientation, score)."""
    m = P.shape[0]
    cols = np.array([A.col(p) for p in range(lo, hi + 1) if p != 0])
    out = [[] for _ in range(S.shape[0])]
    for ori, Q in ((1, P), (-1, Pf.rc(P))) if both else ((1, P),):
        ext = np.vstack([Q.T, np.full((1, m), -99.0)]).astype(np.float32)
        sc = np.zeros((S.shape[0], len(cols)), np.float32)
        for c in range(m):
            sc += ext[:, c][S[:, cols + c]]
        ii, jj = np.nonzero(sc >= thr)
        for i, j in zip(ii, jj):
            out[i].append((int(cols[j]), ori, float(sc[i, j])))
    return out


def iupac_detail(S, pat, lo, hi):
    out = [[] for _ in range(S.shape[0])]
    lut = np.zeros((len(pat), 5), bool)
    for c, ch in enumerate(pat):
        for b in A.IUPAC[ch]:
            lut[c, "ACGT".index(b)] = True
    for p in range(lo, hi + 1):
        if p == 0:
            continue
        j = A.col(p)
        ok = np.ones(S.shape[0], bool)
        for c in range(len(pat)):
            ok &= lut[c][S[:, j + c]]
        for i in np.nonzero(ok)[0]:
            out[i].append((j, 1, None))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--families", default="hpsad_families_v7_down100.tsv")
    ap.add_argument("--tags", default="_f1,_f2,_f3,_f4,_f5,_f6")
    ap.add_argument("--down-tag", default="_f")
    ap.add_argument("--min-reps", type=int, default=3)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    root = Path(a.root)
    mdir = root / "datasets" / "motifs"
    out = Path(a.out) if a.out else root / "datasets" / "profiles"
    src = ["EPDnew", "CAGE"]
    ids_u, seqs_u, w_u, _ = M.load_cohort(root, "C1", "!pseudogene", False, sources=src, families=a.families)
    ids, seqs, w, cls = M.load_cohort(root, "C1", "!pseudogene", False, sources=src, families=a.families,
                                      segments="hpsad_segments_v7_down100.tsv.gz")
    assert ids == ids_u
    Su, S = M.encode(seqs_u), M.encode(seqs)
    print(f"{len(ids):,} promoters")

    # gene and coordinate details, full composite sequence (-2000..+100, repeat-masked bases in lower case)
    info = {}
    csv.field_size_limit(10 ** 8)
    want = set(ids)
    with gzip.open(root / "datasets" / "cohorts" / "hpsad_segments_v7_down100.tsv.gz", "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["anchor_id"] in want and r["cohort"] == "C1":
                info[r["anchor_id"]] = r

    els = Pf.build_elements(mdir, "C1_evidence", a.tags.split(","), a.min_reps)
    hits, names, motif_of = [], [], {}
    taken = set()
    for short, words, lo, hi, strand in A.UP:
        reps = []
        for wd in words:
            reps = [e for k, e in enumerate(els) if k not in taken and M.word_strand(e["rep"]["cons"], [wd])]
            if reps:
                break
        taken |= {k for k, e in enumerate(els) if any(e is r for r in reps)}
        per = [[] for _ in ids]
        for e in reps:
            r = e["rep"]
            flip = M.word_strand(r["cons"], words) < 0
            pwm = r["pwm"][::-1, ::-1] if flip else r["pwm"]
            cons = M.consensus(pwm) if hasattr(M, "consensus") else r["cons"]
            P = np.log2((pwm + 0.01) / 1.04 / 0.25)
            thr = Pf.calibrate(Su, P, r["d"], r["w"], r["n"])
            for i, lst in enumerate(pwm_detail(S, P, thr, lo, hi, strand == "both")):
                per[i] += [(c, o, s, thr, cons, P.shape[0]) for c, o, s in lst]
        hits.append(per); names.append(short)
        print(f"  {short:<7} <- " + (", ".join(e["rep"]["cons"] for e in reps) or "no element"))
    for short, pat, lo, hi in A.CONS:
        hits.append([[(c, o, None, None, pat, len(pat)) for c, o, _ in lst] for lst in iupac_detail(S, pat, lo, hi)])
        names.append(short)
    rows = [r for r in csv.DictReader(open(mdir / f"hpsad_motifs_C1ev_down_zones{a.down_tag}.tsv", encoding="utf-8"),
                                      delimiter="\t")
            if r["pass_fdr"] == "True" and M.word_strand(r["consensus"], A.YY1[1]) == 1
            and A.YY1[2] <= int(r["position"]) <= A.YY1[3] + 10]
    r = max(rows, key=lambda x: float(x["score"]))
    P = np.log2((np.array([[float(x) for x in c.split(",")] for c in r["pwm_ACGT"].split(";")]) + 0.01) / 1.04 / 0.25)
    thr = Pf.calibrate(S, P, int(r["start_depth"]), int(r["window"]), int(r["instances"]))
    hits.append([[(c, o, s, thr, r["consensus"], P.shape[0]) for c, o, s in lst]
                 for lst in pwm_detail(S, P, thr, A.YY1[2], A.YY1[3], False)])
    names.append("YY1")

    # check: presence/absence must equal the published calls
    C = np.array([[len(hits[k][i]) > 0 for k in range(len(names))] for i in range(len(ids))])
    ref = {row["anchor_id"]: row for row in csv.DictReader(open(root / "datasets" / "profiles" / "architecture_calls_C1ev_v7.tsv",
                                                                 encoding="utf-8"), delimiter="\t")}
    mism = sum(int(C[i, k]) != int(ref[ids[i]][nm]) for i in range(len(ids)) for k, nm in enumerate(names))
    print(f"check against architecture_calls_C1ev_v7.tsv: {len(ref)} promoters in file, {mism} mismatching calls")
    if mism or set(ref) != set(ids):
        sys.exit("calls do not reproduce the published table - stopping")

    # tables
    prom_rows, el_rows = [], []
    for i, aid in enumerate(ids):
        g = info[aid]
        chrom, strand, tss = chrom_name(g["chrom"]), g["strand"], int(g["anchor_tss"])
        full = g["sequence_softmasked"]                              # -2000..+100
        up, down = full[:-SHIFT], full[-SHIFT:]
        lo_g, hi_g = sorted((genomic(tss, strand, -len(up)), genomic(tss, strand, SHIFT - 1)))
        present = [nm for k, nm in enumerate(names) if C[i, k]]
        arch = " + ".join(present) if present else "(none)"
        prom_rows.append(dict(
            promoter_id=aid, gene_symbol=g["symbol"], gene_class=cls[i], biotype=g["biotype"],
            chromosome=chrom, strand=strand, start_site_GRCh38=tss, start_site_source=g["anchor_source"],
            region_start_GRCh38=lo_g, region_end_GRCh38=hi_g, family_weight=round(float(w[i]), 4),
            cpg_island_like="yes" if cpg_class(seqs_u[i]) == "cpg" else "no",
            n_elements=len(present), architecture=arch,
            **{f"has_{nm}": int(C[i, k]) for k, nm in enumerate(names) if nm in SUPPORT},
            sequence_core_m300_p100=full[-(300 + SHIFT):], sequence_upstream_m2000_m1=up, sequence_downstream_p1_p100=down,
            ucsc_promoter=ucsc(chrom, lo_g, hi_g),
            ncbi_gene=f"https://www.ncbi.nlm.nih.gov/gene/?term={g['symbol']}%5Bsym%5D+AND+human%5Borgn%5D"))
        for k, nm in enumerate(names):
            for c, o, s, thr, cons, m in sorted(hits[k][i]):
                t0, t1 = offset(c), offset(c + m - 1)
                matched = seqs[i][c:c + m]
                g0, g1 = sorted((genomic(tss, strand, t0), genomic(tss, strand, t1)))
                ctx = seqs[i][max(0, c - 10):c] + "[" + matched + "]" + seqs[i][c + m:c + m + 10]
                el_rows.append(dict(
                    promoter_id=aid, gene_symbol=g["symbol"], element=nm,
                    support=SUPPORT[nm],
                    motif_consensus=cons, match_strand="sense" if o == 1 else "antisense",
                    start_position=rel_pos(t0), end_position=rel_pos(t1),
                    matched_sequence_gene_strand=matched,
                    matched_sequence_motif_orientation=matched if o == 1 else revcomp(matched),
                    matched_sequence_genome_plus_strand=matched if strand == "+" else revcomp(matched),
                    context_gene_strand=ctx, chromosome=chrom, genomic_start_GRCh38=g0, genomic_end_GRCh38=g1,
                    gene_strand=strand, score=round(s, 2) if s is not None else "",
                    threshold=round(thr, 2) if thr is not None else "exact match",
                    ucsc_element=ucsc(chrom, g0, g1, pad=20)))
    out.mkdir(parents=True, exist_ok=True)
    for nm, rows_ in (("S1_Table_promoters.tsv", prom_rows), ("S1_Table_elements.tsv", el_rows)):
        with open(out / nm, "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows_[0]), delimiter="\t")
            wr.writeheader()
            wr.writerows(rows_)
    print(f"{len(prom_rows):,} promoters, {len(el_rows):,} element matches written to {out}")


if __name__ == "__main__":
    main()
