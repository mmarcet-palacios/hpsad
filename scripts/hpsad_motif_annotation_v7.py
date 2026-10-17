#!/usr/bin/env python3
"""
hpsad_motif_annotation_v7.py  -  HPSAD step 5: elements and identities of discovered motifs

Input : datasets/motifs/hpsad_motifs_<cohort>_zones<tag>.tsv   (from hpsad_motifs_v7.py)
        JASPAR 2026 CORE vertebrates non-redundant PFMs (downloaded once to datasets/raw)

1. IDENTITY  Each motif passing FDR is compared with every JASPAR matrix on both
   strands and at every offset with >= 5 overlapping columns. Similarity = sum over
   overlapping columns of the Pearson correlation between the two probability
   columns (the Tomtom "PCC" score). Significance: the same comparison against
   column-shuffled JASPAR matrices (20 shuffles each) gives a null distribution of
   best scores per target length; p = fraction of null scores >= observed and
   E = p x number of JASPAR matrices. Top 5 matches are reported.
2. ELEMENTS  Motifs are linked when their own matrices match each other on either
   strand: >= 5 aligned columns that are informative (>= 1 bit) in both, with
   mean column correlation >= 0.90 over those columns, or when both have the same
   best JASPAR match at E <= 0.05. Linked motifs form one element. Within an element, motif start
   positions closer than 25 nt form one site.

Outputs (datasets/motifs/):
  hpsad_motif_annotation_<cohort><tag>.tsv   one row per motif: element, site, top JASPAR matches
  hpsad_elements_<cohort><tag>.tsv           one row per element: members, sites, best JASPAR match
Usage:
  python scripts\\hpsad_motif_annotation_v7.py --tag _zones1
"""

import argparse
import csv
import math
import io
import ssl
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

JASPAR_TXT = "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
JASPAR_ZIP = "https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.zip"
LOCAL = "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
MIN_OVERLAP = 5
N_SHUFFLE = 20
ELEMENT_PCC = 0.90        # mean column correlation over informative columns
ELEMENT_OVERLAP = 5       # informative aligned columns required
IC_MIN = 1.0              # bits; a column is informative in a motif above this
SAME_HIT_E = 0.05         # motifs whose best JASPAR match is the same matrix at E <= this also join
SITE_GAP = 25


def ssl_context():
    """verified TLS context; uses the certifi bundle when Python's own store is missing."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def fetch_jaspar(raw):
    path = raw / LOCAL
    if path.exists():
        return path.read_text()
    ctx = ssl_context()
    for url in (JASPAR_TXT, JASPAR_ZIP):
        try:
            print(f"Downloading {url}")
            data = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "HPSAD/1.0"}),
                                          timeout=120, context=ctx).read()
            if url.endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(data)) as z:
                    text = "".join(z.read(n).decode() for n in z.namelist() if n.endswith(".jaspar") or n.endswith(".txt"))
            else:
                text = data.decode()
            path.write_text(text)
            return text
        except Exception as e:                       # try the next form
            print(f"  failed: {e}")
    raise SystemExit("Could not download JASPAR; place the file in datasets/raw as " + LOCAL)


def parse_jaspar(text):
    mats, name, rows = [], None, {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name and len(rows) == 4:
                mats.append((name, to_prob(rows)))
            name, rows = line[1:].strip(), {}
        else:
            base = line[0]
            nums = line[line.index("[") + 1:line.index("]")].split() if "[" in line else line[1:].split()
            rows[base] = [float(x) for x in nums]
    if name and len(rows) == 4:
        mats.append((name, to_prob(rows)))
    return mats


def to_prob(rows):
    m = np.array([rows[b] for b in "ACGT"], dtype=float).T          # (len, 4)
    m = m + 0.25
    return m / m.sum(axis=1, keepdims=True)


def rc(p):
    return p[::-1, ::-1]


def col_pcc(a, b):
    """Pearson correlation of every column of a with every column of b: (la, lb)."""
    az = a - a.mean(axis=1, keepdims=True)
    bz = b - b.mean(axis=1, keepdims=True)
    num = az @ bz.T
    den = np.sqrt((az ** 2).sum(axis=1))[:, None] * np.sqrt((bz ** 2).sum(axis=1))[None, :]
    return np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)


def best_alignment(q, t, min_overlap=MIN_OVERLAP):
    """best sum-of-PCC over offsets and strands; returns (score, overlap, offset, strand, mean)."""
    best = (-np.inf, 0, 0, "+", 0.0)
    for strand, tt in (("+", t), ("-", rc(t))):
        C = col_pcc(q, tt)
        lq, lt = C.shape
        for off in range(-(lt - min_overlap), lq - min_overlap + 1):
            ov = min(lq, off + lt) - max(0, off)
            if ov < min_overlap:
                continue
            s = float(np.trace(C, offset=-off))          # sum of C[i, i - off]
            if s > best[0]:
                best = (s, ov, off, strand, s / ov)
    return best


def ic(p):
    return 2 + (p * np.log2(np.maximum(p, 1e-9))).sum(axis=1)


def element_link(a, b):
    """True when a and b align (either strand, any offset) with >= ELEMENT_OVERLAP columns
    informative in both and mean column correlation >= ELEMENT_PCC over those columns."""
    ia = ic(a) >= IC_MIN
    for bb in (b, rc(b)):
        ib = ic(bb) >= IC_MIN
        C = col_pcc(a, bb)
        la, lb = C.shape
        for off in range(-(lb - 1), la):
            idx = [(i, i - off) for i in range(max(0, off), min(la, off + lb)) if ia[i] and ib[i - off]]
            if len(idx) >= ELEMENT_OVERLAP and np.mean([C[i, j] for i, j in idx]) >= ELEMENT_PCC:
                return True
    return False


def parse_pwm(field):
    return np.array([[float(x) for x in row.split(",")] for row in field.split(";")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--cohort", default="C1")
    ap.add_argument("--tag", default="_zones1")
    ap.add_argument("--seed", type=int, default=5)
    a = ap.parse_args()
    root = Path(a.root)
    mdir = root / "datasets" / "motifs"
    rng = np.random.default_rng(a.seed)

    db = parse_jaspar(fetch_jaspar(root / "datasets" / "raw"))
    print(f"JASPAR matrices: {len(db)}")
    rows = [r for r in csv.DictReader(open(mdir / f"hpsad_motifs_{a.cohort}_zones{a.tag}.tsv", encoding="utf-8"),
                                      delimiter="\t") if r["pass_fdr"] == "True"]
    for r in rows:
        r["pwm"] = parse_pwm(r["pwm_ACGT"])
        r["id"] = f"{r['zone']}:M{int(r['rank']):02d}"
    print(f"Motifs passing FDR: {len(rows)}")
    if not rows:
        print("Nothing to annotate.")
        return

    # null: column-shuffled JASPAR matrices
    shuffled = [(n, m[rng.permutation(len(m))]) for n, m in db for _ in range(N_SHUFFLE)]

    lens = {name: len(m) for name, m in db}
    for r in rows:
        q = r["pwm"]
        obs = [(best_alignment(q, m), name) for name, m in db]
        null_by_len = defaultdict(list)
        for (name, m) in shuffled:
            null_by_len[len(m)].append(best_alignment(q, m)[0])
        gum = {}
        for ln, v in null_by_len.items():            # Gumbel fit (method of moments) per target length
            v = np.asarray(v)
            v = v[np.isfinite(v)]                      # -inf = no alignment with >= MIN_OVERLAP columns (target or query too short)
            if len(v) < 10:
                gum[ln] = None
                continue
            beta = max(float(v.std()) * math.sqrt(6) / math.pi, 1e-6)
            gum[ln] = (float(v.mean()) - 0.5772 * beta, beta)
        hits = []
        for (sc, ov, off, strand, mean), name in obs:
            g = gum.get(lens[name])
            if g is None or not np.isfinite(sc):
                p = 1.0                                  # no valid alignment: not a match
            else:
                mu, beta = g
                p = -math.expm1(-math.exp(-(sc - mu) / beta))  # P(best null score >= sc)
            hits.append((p * len(db), sc, mean, ov, strand, name))
        hits.sort(key=lambda x: (x[0], -x[1]))
        r["matches"] = hits[:5]
        top = hits[0]
        print(f"  {r['id']:<12} {r['consensus']:<18} -{r['start_depth']:<5} best {top[5]:<28} "
              f"E={top[0]:.2g} PCC={top[2]:.2f} over {top[3]} ({top[4]})")

    # elements: link motifs whose matrices match each other
    n = len(rows)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            mi, mj = rows[i]["matches"][0], rows[j]["matches"][0]
            same_hit = mi[5] == mj[5] and mi[0] <= SAME_HIT_E and mj[0] <= SAME_HIT_E
            if same_hit or element_link(rows[i]["pwm"], rows[j]["pwm"]):
                parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(rows[i])
    elements = sorted(groups.values(), key=lambda g: -max(float(r["score"]) for r in g))

    out_m, out_e = [], []
    for ei, g in enumerate(elements, 1):
        eid = f"EL{ei:02d}"
        g.sort(key=lambda r: int(r["start_depth"]))
        sites, cur = [], [g[0]]
        for r in g[1:]:
            if int(r["start_depth"]) - int(cur[-1]["start_depth"]) < SITE_GAP and r["zone"] == cur[-1]["zone"]:
                cur.append(r)
            else:
                sites.append(cur)
                cur = [r]
        sites.append(cur)
        best = min((m for r in g for m in r["matches"]), key=lambda x: x[0])
        votes = defaultdict(int)
        for r in g:
            if r["matches"][0][0] <= 0.05:
                votes[r["matches"][0][5]] += 1
        for si, s in enumerate(sites, 1):
            for r in s:
                r["element"], r["site"] = eid, f"{eid}.{si}"
        rep = max(g, key=lambda r: float(r["score"]))
        out_e.append({"element": eid, "n_motifs": len(g), "representative": rep["consensus"],
                      "best_jaspar": best[5], "best_E": f"{best[0]:.2g}", "best_PCC": round(best[2], 2),
                      "top_hit_votes": "; ".join(f"{k} ({v})" for k, v in sorted(votes.items(), key=lambda x: -x[1])),
                      "sites": "; ".join(f"-{min(int(r['start_depth']) for r in s)}..-{max(int(r['start_depth']) for r in s)}"
                                         f" ({len(s)})" for s in sites),
                      "zones": ",".join(sorted({r["zone"] for r in g})),
                      "members": ", ".join(f"{r['id']} {r['consensus']} -{r['start_depth']}" for r in g),
                      "max_score": max(float(r["score"]) for r in g),
                      "max_ratio": max(float(r["ratio"]) for r in g)})
    for r in rows:
        out_m.append({"motif": r["id"], "element": r["element"], "site": r["site"], "consensus": r["consensus"],
                      "start_depth": r["start_depth"], "window": r["window"], "ratio": r["ratio"], "score": r["score"],
                      "instances": r["instances"], "known_rule": r["known"],
                      **{f"jaspar{i + 1}": f"{m[5]} | E={m[0]:.2g} | PCC={m[2]:.2f} | ov={m[3]} | {m[4]}"
                         for i, m in enumerate(r["matches"])}})

    def write(path, rs):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rs[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rs)
    write(mdir / f"hpsad_motif_annotation_{a.cohort}{a.tag}.tsv", sorted(out_m, key=lambda x: x["motif"]))
    write(mdir / f"hpsad_elements_{a.cohort}{a.tag}.tsv", out_e)
    print(f"\n{len(rows)} motifs -> {len(elements)} elements")
    for e in out_e:
        print(f"  {e['element']} {e['representative']:<18} {e['best_jaspar']:<28} E={e['best_E']:<8} "
              f"sites {e['sites']}  [{e['n_motifs']} motifs]")


if __name__ == "__main__":
    main()
