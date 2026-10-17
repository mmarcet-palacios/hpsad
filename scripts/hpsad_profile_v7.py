#!/usr/bin/env python3
"""
hpsad_profile_v7.py  -  where along the promoter does each reproducible element sit?

Discovery finds an element once, at one window, and then masks it; a second site of the same
element (e.g. a second CCAAT box) has to compete with the first one sitting in its flank. Here
each element is instead scanned along every promoter at every position, on both strands, and
compared with the SAME scan on the ten shuffled cohorts of the run (dinucleotide shuffles in
50-nt blocks, which keep local composition). The expectation at each position is the mean of
the ten shuffles - no flanks, no masking, no random draws.

Elements   motifs passing FDR in the replicate runs, grouped by PWM shape alone (any position;
           mean column Pearson r >= 0.85 over >= 5 informative columns, offsets -6..+6), with
           reverse-complement pairs merged (e.g. CCAAT and ATTGG). Kept when the element passes
           in >= --min-reps replicates. Copy artefacts (O/E > 20 and < 60 carriers) are skipped.
           The representative PWM is the highest-scoring passing motif.
Threshold  per element, calibrated like hpsad_coverage_v7.py: in its own discovery window
           (sense strand), the number of promoters scoring >= t equals its instance count.
Profile    O(d) = family-weighted number of promoters with a match starting at upstream position d
           (either strand; a start matching both strands counts once); E(d) = mean of the ten
           shuffled cohorts. Sense and antisense profiles are kept separately.
Summits    O and E are summed over 11 nt around each position (E from the shuffle density averaged
           over 51 nt, so rare matrices are not divided by near-zero counts); summits are local maxima
           of log2((O+0.5)/(E+0.5)) with a prominence >= 0.4 (a 1.3-fold dip to the next summit), inside
           --zone, scoring above the cut (score = -log10 P(Poisson(E) >= O)). Cut = the highest summit
           score over all elements in the null runs, where each shuffled cohort in turn plays "real"
           against the mean of the other nine.

Usage (from C:\\Promoter\\HPSAD_v7):
    python scripts\\hpsad_profile_v7.py --label C1_evidence --tags _v5s1,_v5s2,_v5s3 --sources EPDnew,CAGE
Outputs (datasets\\profiles):
    profile_<label><tag1>.tsv          element x position: O and E, sense and antisense
    profile_peaks_<label><tag1>.tsv    peaks per element (primary, secondary, ...)
    profile_<label><tag1>_summary.json elements, thresholds, null peak scores, peaks
"""
import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpsad_motifs_v7 as M                                   # noqa: E402
import hpsad_stability_v7 as St                               # noqa: E402

L = M.L
PEAK_WIDTHS = (1, 5, 11, 21, 31)


def rc(p):
    return p[::-1, ::-1]


def similar_any_pos(a, b):
    return St.similar(dict(a, d=0, zone="x"), dict(b, d=0, zone="x"), 10, 0.85)


def build_elements(mdir, label, tags, min_reps):
    motifs = []
    for t in tags:
        ann = {r["motif"]: r for r in csv.DictReader(open(mdir / f"hpsad_motif_annotation_{label}{t}.tsv", encoding="utf-8"),
                                                     delimiter="\t")}
        win = {(r["zone"], int(r["rank"])): int(r["window"])
               for r in csv.DictReader(open(mdir / f"hpsad_motifs_{label}_zones{t}.tsv", encoding="utf-8"), delimiter="\t")}
        for m in St.load(mdir / f"hpsad_motifs_{label}_zones{t}.tsv"):
            m["w"] = win[(m["zone"], m["rank"])]
            if not m["passed"] or (m["oe"] > 20 and m["n"] < 60):
                continue
            mid = f"{m['zone']}:M{m['rank']:02d}"
            j = ann.get(mid, {}).get("jaspar1", "")
            m.update(tag=t, mid=mid, jaspar=j.split("|")[0].split("\t")[-1].strip() if j else "",
                     pcc=j.split("PCC=")[1].split()[0] if "PCC=" in j else "")
            motifs.append(m)
    els = []                                        # each: list of (motif, orientation +1/-1)
    for m in sorted(motifs, key=lambda x: -x["score"]):
        for e in els:
            rep = e[0][0]
            if similar_any_pos(m, rep):
                e.append((m, 1))
                break
            if similar_any_pos(dict(m, pwm=rc(m["pwm"])), rep):
                e.append((m, -1))
                break
        else:
            els.append([(m, 1)])
    out = []
    for e in els:
        reps_passing = sorted({m["tag"] for m, _ in e})
        if len(reps_passing) < min_reps:
            continue
        rep = e[0][0]
        name = rep["known"] or rep["jaspar"] or rep["cons"]
        out.append(dict(name=name, rep=rep, members=e, n_reps=len(reps_passing),
                        positions=sorted(m["d"] for m, _ in e)))
    return out


def scan_hits(S, logodds, thr):
    """bool (n, L-m+1): window starting at column j scores >= thr (column j <-> upstream position L-j)."""
    m = logodds.shape[0]
    ext = np.vstack([logodds.T, np.full((1, m), -99.0)]).astype(np.float32)      # (5, m); code 4 = masked
    J = L - m + 1
    sc = np.zeros((S.shape[0], J), np.float32)
    for c in range(m):
        sc += ext[:, c][S[:, c:c + J]]
    return sc >= thr


def calibrate(S, logodds, d0, w, inst):
    m = logodds.shape[0]
    ext = np.vstack([logodds.T, np.full((1, m), -99.0)])
    lo, hi = max(m, d0 - w // 2 - 3), min(L, d0 + w // 2 + 3)
    best = np.full(S.shape[0], -1e9)
    for d in range(lo, hi + 1):
        j = L - d
        best = np.maximum(best, ext[S[:, j:j + m], np.arange(m)[None, :]].sum(1))
    return float(np.sort(best)[::-1][min(inst, len(best)) - 1])


def profiles(S, w, P, thr):
    hs = scan_hits(S, P, thr)
    ha = scan_hits(S, rc(P), thr)
    to_depth = lambda v: v[::-1]                   # column j -> index of depth d = L - j (index 0 = depth m ... )
    return (to_depth(w @ (hs | ha)), to_depth(w @ hs), to_depth(w @ ha),
            to_depth((hs | ha).sum(0).astype(float)))


SMOOTH = 11            # summit window (nt)
E_SPAN = 51            # expectation: shuffle density averaged over this span (shuffles carry no signal below ~50 nt)
MIN_PROM = 0.4         # log2 prominence of a summit (1.3-fold dip to the neighbouring summit)


def running_sum(x, w):
    c = np.concatenate([[0.0], np.cumsum(x)])
    h = w // 2
    i = np.arange(len(x))
    return c[np.clip(i + h + 1, 0, len(x))] - c[np.clip(i - h, 0, len(x))]


def curves(O, E):
    """O and E summed over the SMOOTH-nt window centred on each position; E from the E_SPAN-nt density."""
    Ow = running_sum(O, SMOOTH)
    dens = running_sum(E, E_SPAN) / np.maximum(running_sum(np.ones_like(E), E_SPAN), 1)
    Ew = dens * running_sum(np.ones_like(E), SMOOTH)
    return Ow, Ew, np.log2((Ow + 0.5) / (Ew + 0.5))


def summits(O, E, depths, zone, cut, max_peaks=6):
    from scipy.signal import find_peaks
    Ow, Ew, lr = curves(O, E)
    sc = M.poisson_score(Ow, np.maximum(Ew, 1e-9))
    idx, prop = find_peaks(lr, prominence=MIN_PROM)
    keep = [(i, p) for i, p in zip(idx, prop["prominences"])
            if zone[0] <= depths[i] <= zone[1] and sc[i] >= cut and Ow[i] > Ew[i]]
    keep.sort(key=lambda x: -sc[x[0]])
    return [dict(score=float(sc[i]), centre=int(depths[i]), width=SMOOTH, O=float(Ow[i]), E=float(Ew[i]),
                 log2_prominence=round(float(p), 2)) for i, p in keep[:max_peaks]]


def best_null_score(O, E, depths, zone):
    """highest summit score in a null profile (no prominence or score requirement beyond being a summit)."""
    from scipy.signal import find_peaks
    Ow, Ew, lr = curves(O, E)
    sc = M.poisson_score(Ow, np.maximum(Ew, 1e-9))
    idx, _ = find_peaks(lr, prominence=MIN_PROM)
    idx = [i for i in idx if zone[0] <= depths[i] <= zone[1] and Ow[i] > Ew[i]]
    return float(max(sc[idx])) if idx else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--label", default="C1_evidence")
    ap.add_argument("--tags", default="_v5s1,_v5s2,_v5s3")
    ap.add_argument("--cohort", default="C1")
    ap.add_argument("--gene-class", default="!pseudogene")
    ap.add_argument("--sources", default="EPDnew,CAGE")
    ap.add_argument("--min-reps", type=int, default=2, help="element must pass in at least this many replicates")
    ap.add_argument("--zone", default="1,300", help="peak search range (upstream positions)")
    ap.add_argument("--nulls", type=int, default=10)
    ap.add_argument("--families", default="hpsad_families_v7_down100.tsv", help="family weights file in datasets/families")
    ap.add_argument("--suffix", default="", help="appended to output names (e.g. _fam2)")
    a = ap.parse_args()
    root = Path(a.root)
    mdir, out = root / "datasets" / "motifs", root / "datasets" / "profiles"
    out.mkdir(parents=True, exist_ok=True)
    tags = a.tags.split(",")
    zone = tuple(int(x) for x in a.zone.split(","))

    ids, seqs, w, cls = M.load_cohort(root, a.cohort, a.gene_class, False,
                                      sources=a.sources.split(",") if a.sources else None, families=a.families)
    S = M.encode(seqs)
    shuf = []
    for r in range(a.nulls):
        p = mdir / "nulls" / f"null50_{a.label}{tags[0]}_rep{r}.txt.gz"
        ns = gzip.open(p, "rt").read().split("\n")
        if len(ns) != len(seqs):
            sys.exit(f"{p.name}: {len(ns)} sequences, cohort has {len(seqs)} - wrong cohort options?")
        shuf.append(M.encode(ns))
    print(f"{a.label}: {len(seqs):,} promoters ({w.sum():,.1f} independent); {len(shuf)} shuffled cohorts from {tags[0]}")

    els = build_elements(mdir, a.label, tags, a.min_reps)
    print(f"{len(els)} elements passing in >= {a.min_reps} of {len(tags)} replicates")
    prof_rows, results = [], []
    null_first = []
    for k, el in enumerate(els, 1):
        rep = el["rep"]
        P = np.log2((rep["pwm"] + 0.01) / 1.04 / 0.25)
        m = P.shape[0]
        thr = calibrate(S, P, rep["d"], rep["w"], rep["n"])
        depths = np.arange(m, L + 1)                                # index i <-> upstream position m + i
        Oa, Os, Oas, nO = profiles(S, w, P, thr)
        sh = [profiles(Z, w, P, thr) for Z in shuf]
        Ea = np.mean([x[0] for x in sh], axis=0)
        Es = np.mean([x[1] for x in sh], axis=0)
        Eas = np.mean([x[2] for x in sh], axis=0)
        for r in range(len(sh)):                                    # null: shuffle r against the other nine
            En = np.mean([x[0] for q, x in enumerate(sh) if q != r], axis=0)
            null_first.append(best_null_score(sh[r][0], En, depths, zone))
        el.update(P=P, thr=thr, depths=depths, Oa=Oa, Os=Os, Oas=Oas, nO=nO, Ea=Ea, Es=Es, Eas=Eas)
        print(f"  [{k}/{len(els)}] {el['name']:<22} {rep['cons']:<18} thr {thr:.2f}  "
              f"passes in {el['n_reps']} reps at {el['positions']}")
    cut = max(null_first) if null_first else 0.0
    print(f"\nsummit cut = highest summit score in {len(null_first)} null profiles: {cut:.1f}")

    peak_rows = []
    for el in els:
        rep, depths = el["rep"], el["depths"]
        peaks = summits(el["Oa"], el["Ea"], depths, zone, cut)
        for i, p in enumerate(peaks, 1):
            h = p["width"] // 2
            sel = (depths >= p["centre"] - h) & (depths <= p["centre"] + h)
            os_, oa_ = el["Os"][sel].sum(), el["Oas"][sel].sum()
            es_, ea_ = el["Es"][sel].sum(), el["Eas"][sel].sum()
            peak_rows.append(dict(element=el["name"], consensus=rep["cons"], peak=i, centre=-p["centre"],
                                  start=-(p["centre"] + h), end=-(p["centre"] - h), width=p["width"],
                                  O=round(p["O"], 1), E=round(p["E"], 1), ratio=round(p["O"] / max(p["E"], 1e-9), 2),
                                  score=round(p["score"], 1), log2_prominence=p["log2_prominence"],
                                  O_weighted=round(p["O"], 1), sense_fraction=round(os_ / max(os_ + oa_, 1e-9), 2),
                                  sense_fraction_expected=round(es_ / max(es_ + ea_, 1e-9), 2),
                                  promoters=int(el["nO"][sel].sum())))
        results.append(dict(name=el["name"], consensus=rep["cons"], jaspar=rep["jaspar"], pcc=rep["pcc"],
                            n_reps=el["n_reps"], discovery_positions=el["positions"], threshold=round(el["thr"], 3),
                            representative=f"{rep['tag']} {rep['mid']}", peaks=[r for r in peak_rows if r["element"] == el["name"]]))
        for i, d in enumerate(depths):
            if d <= 1000:
                prof_rows.append(dict(element=el["name"], position=-int(d),
                                      O=round(float(el["Oa"][i]), 2), E=round(float(el["Ea"][i]), 2),
                                      O_sense=round(float(el["Os"][i]), 2), O_antisense=round(float(el["Oas"][i]), 2),
                                      E_sense=round(float(el["Es"][i]), 2), E_antisense=round(float(el["Eas"][i]), 2)))

    stem = f"{a.label}{tags[0]}{a.suffix}"
    with open(out / f"profile_peaks_{stem}.tsv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(peak_rows[0]) if peak_rows else ["element"], delimiter="\t")
        wr.writeheader()
        wr.writerows(peak_rows)
    with open(out / f"profile_{stem}.tsv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(prof_rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(prof_rows)
    json.dump({"label": a.label, "replicates": tags, "promoters": len(seqs), "independent": round(float(w.sum()), 1),
               "zone": zone, "peak_cut": cut, "null_first_peak_scores": null_first, "elements": results},
              open(out / f"profile_{stem}_summary.json", "w"), indent=1, default=float)
    print(f"\n{'element':<16}{'consensus':<18}{'#':>2}{'centre':>7}{'O/E':>7}{'score':>7}{'dip':>5}{'O_w':>7}{'prom':>6}{'sense':>7}{'exp':>6}")
    for r in peak_rows:
        print(f"{r['element'][:15]:<16}{r['consensus']:<18}{r['peak']:>2}{r['centre']:>7}{r['ratio']:>7}{r['score']:>7}"
              f"{r['log2_prominence']:>5}{r['O_weighted']:>7}{r['promoters']:>6}{r['sense_fraction']:>7}{r['sense_fraction_expected']:>6}")
    print(f"\nwritten {out / f'profile_peaks_{stem}.tsv'} and profile_{stem}.tsv")


if __name__ == "__main__":
    main()
