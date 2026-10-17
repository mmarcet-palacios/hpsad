import json, math, csv, os, re
D = json.load(open("data.json"))
ROOT = "out/project"
os.makedirs(ROOT, exist_ok=True)
INK, MUTED, LINE, BG, PANEL = "#23231f", "#6b675c", "#d9d4c7", "#f7f5ef", "#ffffff"
BLUE, ORANGE, RED, GREEN, GREY = "#3f6fa8", "#d9731f", "#b5452f", "#3d8a5f", "#a9a497"
FONTS = ('<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&amp;'
         'family=IBM+Plex+Mono&amp;family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&amp;display=swap" rel="stylesheet">')
W = 1040
SVGW = 928
LEG = json.load(open("legends.json"))
CUR = None          # manuscript figure number of the board being built (legend shown in the footer box)


def legend_px():
    if not CUR:
        return 0
    t, b = LEG[CUR]
    return int(math.ceil((len(t) + len(b)) / 122) * 19.5 + 44)


HFIX = json.load(open("heights_fix.json")) if os.path.exists("heights_fix.json") else {}


def page(title, h, body):
    h = h + legend_px()
    if CUR and f"Fig{CUR}.dc.html" in HFIX:
        h = HFIX[f"Fig{CUR}.dc.html"] + 60     # measured natural height of the content + margin
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="./support.js"></script>
</head>
<body>
<x-dc>
<helmet>
{FONTS}
<style>body{{margin:0;font-family:'IBM Plex Sans',sans-serif;color:{INK};background:{BG}}}</style>
</helmet>
<div style="width: {W}px; height: {h}px; box-sizing: border-box; padding: 36px 56px; display: flex; flex-direction: column; gap: 14px; background: {BG}">
{body}
</div>
</x-dc>
<script type="text/x-dc" data-dc-script data-props='{{"$preview":{{"width":{W},"height":{h}}}}}'>
class Component extends DCLogic {{
  renderVals() {{ return {{}}; }}
}}
</script>
</body>
</html>
"""


def header(num, title, message, status, status_col="#e6efe0", status_ink="#2f5d3a"):
    return f"""<div style="display: flex; align-items: center; gap: 12px">
<div style="font-size: 12px; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Figure {num}</div>
<div style="font-size: 12px; padding: 3px 10px; border-radius: 999px; background: {status_col}; color: {status_ink}">{status}</div>
</div>
<div style="font-family: 'Source Serif 4', serif; font-size: 25px; font-weight: 600; line-height: 1.25">{title}</div>
<div style="font-size: 15px; line-height: 1.45; color: #3b3a34; max-width: 900px">{message}</div>"""


def footer(sources, notes):
    n = "".join(f"<div>· {x}</div>" for x in notes)
    leg = ""
    if CUR:
        t, b = LEG[CUR]
        leg = (f'<div style="grid-column: 1 / -1; padding-bottom: 12px; border-bottom: 1px solid #d6d0bf">'
               f'<div style="font-weight: 600; margin-bottom: 4px">Legend (manuscript draft)</div>'
               f'<div style="color: #3b3a34"><b>{t}</b> {b}</div></div>')
    return f"""<div style="margin-top: auto; display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 14px 24px; padding: 14px 18px; background: #ece8dc; border-radius: 10px; font-size: 13px; line-height: 1.5">
{leg}
<div><div style="font-weight: 600; margin-bottom: 4px">Data</div><div style="font-family: 'IBM Plex Mono', monospace; font-size: 12px; color: #3b3a34">{sources}</div></div>
<div><div style="font-weight: 600; margin-bottom: 4px">To decide</div>{n}</div>
</div>"""


def svg(h, inner, label):
    return (f'<svg width="{SVGW}" height="{h}" viewBox="0 0 {SVGW} {h}" role="img" aria-label="{label}" '
            f'style="background: {PANEL}; border: 1px solid {LINE}; border-radius: 8px">'
            f'<g font-family="IBM Plex Sans, sans-serif" font-size="12" fill="{INK}">{inner}</g></svg>')


def T(x, y, s, size=12, fill=INK, anchor="start", weight=None, style=None, family=None):
    a = f' font-weight="{weight}"' if weight else ""
    a += f' font-style="{style}"' if style else ""
    a += f' font-family="{family}"' if family else ""
    s = re.sub(r"<(?=[\s\d=])", "&lt;", str(s))
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}"{a}>{s}</text>'


def R(x, y, w, h, fill, stroke="none", sw=1, rx=0, extra=""):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w,0):.1f}" height="{max(h,0):.1f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}"{extra}></rect>'


def Lн(x1, y1, x2, y2, stroke=LINE, sw=1, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{sw}"{d}></line>'


L = Lн


def panel_label(x, y, letter, text):
    return T(x, y, letter, 18, weight="600") + T(x + 24, y, text, 12, MUTED)


def axis_y(x, y0, y1, vmax, ticks, fmt=lambda v: f"{v:g}", title=None, log=False, vmin=None):
    out = L(x, y0, x, y1, "#8a8578")
    for v in ticks:
        if log:
            yy = y0 - (math.log10(v) - math.log10(vmin)) / (math.log10(vmax) - math.log10(vmin)) * (y0 - y1)
        else:
            yy = y0 - v / vmax * (y0 - y1)
        out += L(x - 4, yy, x, yy, "#8a8578") + T(x - 7, yy + 4, fmt(v), 10, MUTED, "end")
    if title:
        out += f'<text x="{x-38:.1f}" y="{(y0+y1)/2:.1f}" font-size="10.5" fill="{MUTED}" text-anchor="middle" transform="rotate(-90 {x-38:.1f} {(y0+y1)/2:.1f})">{title}</text>'
    return out


def pct(v):
    return f"{v*100:.0f}%" if v >= 0.1 else f"{v*100:.1f}%"


# ───────────────────────── Main ─────────────────────────
def main_board():
    arc = [
        ("Goal", "Map positioned sequence elements across all human promoters, de novo, with each promoter counted once."),
        ("Problem 1: the anchor", "Gene 5′ edges miss the experimental start site by ≥ 100 nt for 31% of genes (Fig 1)."),
        ("Problem 2: shared ancestry", "Duplicated families look like perfect positioned elements; detect and count them once (Fig 2)."),
        ("A calibrated statistic", "Each word against its own flanks; recovers TATA, GC and CCAAT boxes and 60/60 planted motifs (Fig 3)."),
        ("Know the limits", "Motifs carried by ~140 promoters at a fixed position are found; fuzzy ones need ~280+ (Fig 4)."),
        ("What promoters hold", "Promoters with an EPDnew or CAGE start site: 7 elements pass in all 6 replicates and 3 more in at least 3; their positional profiles show the layout — TATA at −31, a crowded slot at −53 to −55, a second slot at the TATA position in TATA-less promoters, CCAAT spread over −34 to −160 (Fig 5)."),
        ("The start site", "Extending to +100 adds the Initiator, TCT (ribosomal-protein genes), a weak DPE at +28 in TATA-less promoters and YY1 at +10 to +21 (Fig 6)."),
        ("Shared across cohorts", "Same grammar, same positions in truncated, head-to-head and nested promoters; a few architectures with a TATA pole and a CpG pole (Fig 7)."),
        ("Pseudogenes", "The −17 motif is the LINE-1 insertion site, confirmed by target-site duplications (Fig 8)."),
    ]
    a = "".join(f'<div><b>{i+1} · {t}.</b> {s}</div>' for i, (t, s) in enumerate(arc))
    figs = [("Fig 1", "Gene edges vs experimental start sites"), ("Fig 2", "Cohorts, families, independent promoters"),
            ("Fig 3", "Test bench for the positional statistic"), ("Fig 4", "Detection limits of motif discovery"),
            ("Fig 5", "Reproducible elements and their positional profiles"), ("Fig 6", "The start-site region: Initiator, TCT, DPE, YY1"), ("Fig 7", "The grammar across cohorts and gene classes"),
            ("Fig 8", "The −17 motif and target-site duplications")]
    f = "".join(f'<div><b>{n}</b> {t}</div>' for n, t in figs)
    body = f"""<div style="display: flex; flex-direction: column; gap: 6px">
<div style="font-size: 13px; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">HPSAD v7 · Research article · figure storyboard, draft 3 · 14 Oct 2026</div>
<h1 style="margin: 0; font-family: 'Source Serif 4', serif; font-size: 32px; font-weight: 600; line-height: 1.22">[WORKING TITLE] De novo positional discovery across human promoters recovers the core grammar and the footprints of retrotransposition</h1>
</div>
<div style="display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 28px">
<div style="display: flex; flex-direction: column; gap: 10px">
<div style="font-family: 'Source Serif 4', serif; font-size: 20px; font-weight: 600">Story arc</div>
<div style="display: flex; flex-direction: column; gap: 7px; font-size: 14.5px; line-height: 1.42">{a}</div>
</div>
<div style="display: flex; flex-direction: column; gap: 10px">
<div style="font-family: 'Source Serif 4', serif; font-size: 20px; font-weight: 600">Figure list</div>
<div style="display: flex; flex-direction: column; gap: 6px; font-size: 15px; line-height: 1.4">{f}</div>
<div style="font-size: 13px; color: {MUTED}; line-height: 1.45; padding-top: 6px">Supporting (proposed): S1 Table all motifs with JASPAR matches · S2 Table motifs per cohort and gene class · S1 Fig largest families (OR7E, GOLGA, USP17L, SNAR) · S1 File code</div>
<div style="display: flex; gap: 8px; flex-wrap: wrap; padding-top: 4px; font-size: 12px">
<span style="padding: 3px 10px; border-radius: 999px; background: #e6efe0; color: #2f5d3a">Data ready · to draw: Figs 1–8</span>
<span style="padding: 3px 10px; border-radius: 999px; background: #e6efe0; color: #2f5d3a">No-pseudogene runs in (Fig 6A)</span>
</div>
</div>
</div>
<div style="margin-top: auto; padding: 16px 20px; background: #ece8dc; border-radius: 10px; font-size: 14px; line-height: 1.5">
<b>To decide together:</b> Figs 3 and 4 could merge (bench + limits) with one moving to supplementary. · Fig 5 rebuilt 14 Oct around positional profiles (motif-level matrix → Supplementary Fig S2). · Decided 12 Oct: Fig 5 uses promoters with EPDnew or CAGE support (C1 evidence, no pseudogenes); pseudogenes appear in Fig 7 only. · Fig 7: keep B (word enrichment) or move it to supplementary now that C–E show the duplications directly?
<div style="margin-top: 8px"><b>Reminders:</b> (1) Remove round-1 sections from the manuscript. (2) References to add: LINE-1 endonuclease target site and target-site duplications; THAP11/ZNF143 at bidirectional promoters. (3) Figs 7–8 still use the v2/v3 runs; Fig 5 uses v4–v5, Fig 6 the v7 downstream run (15 Oct).</div>
</div>"""
    return page("HPSAD storyboard overview", 860, body)


# ───────────────────────── Fig 1 ─────────────────────────
def fig1():
    rows = list(csv.DictReader(open("/mnt/user-data/uploads/Promoter/HPSAD_v7/datasets/qc/tss_anchor_comparison.csv")))
    rows = [r for r in rows if r["is_main_promoter"] == "True" and r["offset_gene_edge"] != ""]
    out = ""
    for k, (col, lab, key) in enumerate([(GREY, "A", "offset_gene_edge"), (BLUE, "B", "offset_nearest_tx")]):
        x0, y0, pw, ph = 70 + k * 450, 325, 380, 230
        bins = [0] * 100
        outside = 0
        for r in rows:
            v = int(float(r[key]))
            if -500 <= v < 500:
                bins[(v + 500) // 10] += 1
            else:
                outside += 1
        cap = 600
        out += panel_label(x0 - 50, 30, lab, "offset from " + ("gene 5′ edge" if k == 0 else "nearest transcript 5′ end") + " to EPDnew start")
        out += axis_y(x0, y0, y0 - ph, cap, [0, 200, 400, 600], title="genes" if k == 0 else None)
        for i, c in enumerate(bins):
            hh = min(c, cap) / cap * ph
            out += R(x0 + i * pw / 100, y0 - hh, pw / 100 - 0.4, hh, col)
        # bins above the cap: break mark on the bar, labels stacked to the right with leader lines
        capped = [(i, c) for i, c in enumerate(bins) if c > cap]
        for j, (i, c) in enumerate(capped):
            bx = x0 + (i + 0.5) * pw / 100
            out += L(bx - 4, y0 - ph + 8, bx + 4, y0 - ph + 4, INK, 1.6) + L(bx - 4, y0 - ph + 13, bx + 4, y0 - ph + 9, INK, 1.6)
            lx, ly = x0 + pw * 0.62, y0 - ph - 4 - (len(capped) - 1 - j) * 15
            v = i * 10 - 500
            out += L(bx + 1.5, y0 - ph + 1, lx - 4, ly - 4, "#8a8578", 0.8)
            out += T(lx, ly, f"{c:,} genes at {v:+d} to {v + 9:+d}", 10, MUTED)
        out += L(x0, y0, x0 + pw, y0, "#8a8578")
        for v in (-500, -250, 0, 250, 500):
            xx = x0 + (v + 500) / 1000 * pw
            out += L(xx, y0, xx, y0 + 4, "#8a8578") + T(xx, y0 + 16, f"{v:+d}" if v else "0", 10, MUTED, "middle")
        out += T(x0 + pw / 2, y0 + 34, "offset (nt); + = EPDnew start downstream", 10.5, MUTED, "middle")
        out += T(x0 + pw, y0 + 50, f"{outside:,} genes outside ±500 not drawn", 10, MUTED, "end")
    body = header(1, "Fig 1 · Annotated gene edges misplace the start site",
                  "Gene 5′ edges miss the experimental start site by 100 nt or more for nearly a third of genes; the nearest annotated transcript 5′ end is much closer, and EPDnew is used as the anchor.",
                  "Data ready · to draw")
    body += svg(395, out, "Histograms of anchor offsets")
    body += footer("datasets/qc/tss_anchor_comparison.csv<br>tss_anchor_comparison_summary.json<br>(n = 16,442 genes, main EPDnew promoter)",
                   ["Linear y with the 0 bin capped (as sketched) or log y?", "Add a third panel: offset by gene class? EPDnew covers mainly protein-coding genes."])
    return page("Fig 1 sketch", 750, body)


# ───────────────────────── Fig 2 ─────────────────────────
def fig2():
    out = ""
    co = D["fig2_cohorts"]
    cohorts = ["C1", "C2", "C2b", "C4", "C0"]
    cls = [("protein_coding", BLUE, "protein-coding"), ("non_coding", ORANGE, "non-coding"), ("pseudogene", GREY, "pseudogene")]
    x0, y0, ph = 70, 290, 220
    vmax = 30000
    out += panel_label(20, 30, "A", "primary anchors per cohort, by gene class")
    out += axis_y(x0, y0, y0 - ph, vmax, [0, 10000, 20000, 30000], lambda v: f"{v/1000:.0f}k", "anchors")
    for i, c in enumerate(cohorts):
        xx, acc = x0 + 20 + i * 66, 0
        for k, col, _ in cls:
            n = co.get(f"{c}|{k}", 0)
            hh = n / vmax * ph
            out += R(xx, y0 - acc - hh, 40, hh, col)
            acc += hh
        tot = sum(co.get(f"{c}|{k}", 0) for k, _, _ in cls)
        out += T(xx + 20, y0 + 16, c, 11, INK, "middle") + T(xx + 20, y0 - acc - 5, f"{tot:,}", 9.5, MUTED, "middle")
    out += L(x0, y0, x0 + 360, y0, "#8a8578")
    ly = 60
    for k, col, lab in cls:
        out += R(330, ly - 9, 10, 10, col) + T(345, ly, lab, 10.5, MUTED)
        ly += 16
    out += T(330, ly + 4, f"+ C3 alternative {D['fig2_all']['C3']:,}", 10, MUTED) + T(330, ly + 18, f"C3x internal {D['fig2_all']['C3x']:,}", 10, MUTED)
    # B calibration
    bx, by, bw, bh = 520, 290, 340, 220
    out += panel_label(470, 30, "B", "planted copies linked to their source")
    out += axis_y(bx, by, by - bh, 1.0, [0, 0.5, 1.0], lambda v: f"{v:.1f}", "fraction linked")
    cal = D["fig2_cal"]
    ids = ["1.0", "0.95", "0.9", "0.85", "0.82", "0.8", "0.75", "0.7"]
    pts = []
    for i, k in enumerate(ids):
        xx = bx + 14 + (1.0 - float(k)) / 0.3 * 230
        yy = by - cal[k] * bh
        pts.append((xx, yy))
        out += L(xx, by, xx, by + 4, "#8a8578") + T(xx, by + 16, f"{float(k)*100:.0f}", 9.5, MUTED, "middle")
    out += '<polyline points="' + " ".join(f"{x:.1f},{y:.1f}" for x, y in pts) + f'" fill="none" stroke="{BLUE}" stroke-width="2"></polyline>'
    out += "".join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.5" fill="{BLUE}"></circle>' for x, y in pts)
    for j, (k, lab) in enumerate([("shuffled", "shuffled"), ("shared_300nt_chunk", "300-nt chunk")]):
        xx = bx + 268 + j * 40
        out += R(xx, by - 2, 24, 2, ORANGE) + T(xx + 12, by - 8, "0", 10, ORANGE, "middle")
        out += T(xx + 12, by + 16, lab.split()[0], 9, MUTED, "middle")
    out += L(bx, by, bx + bw, by, "#8a8578") + T(bx + 130, by + 34, "true identity (%)", 10.5, MUTED, "middle") + T(bx + 300, by + 34, "controls", 10.5, MUTED, "middle")
    # C families by size bin: share of families by gene class (colour) and arrangement (hatched = dispersed)
    FB = json.load(open("fam_bins.json"))
    binsC = ["2", "3", "4–5", "6–10", "11–20", "21–60"]
    cx, cy, cw, ch = 70, 610, 360, 200
    out += panel_label(20, 360, "C", "families by size: gene class (colour), dispersed copies (hatched)")
    out += ('<defs><pattern id="hatchC" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
            '<line x1="0" y1="0" x2="0" y2="5" stroke="#ffffff" stroke-width="2"></line></pattern></defs>')
    out += axis_y(cx, cy, cy - ch, 1.0, [0, 0.5, 1.0], lambda v: f"{v * 100:.0f}%", "families")
    bwid = (cw - 20) / len(binsC)
    for i, bn in enumerate(binsC):
        c = FB.get(bn, {})
        tot = sum(c.values())
        xx, acc = cx + 10 + i * bwid + 6, 0
        for k, col, _ in cls:
            for arr in ("clustered", "dispersed"):
                n = c.get(f"{k}|{arr}", 0)
                hh = n / tot * ch if tot else 0
                out += R(xx, cy - acc - hh, bwid - 12, hh, col)
                if arr == "dispersed" and hh > 0:
                    out += R(xx, cy - acc - hh, bwid - 12, hh, "url(#hatchC)")
                acc += hh
        out += T(xx + (bwid - 12) / 2, cy + 16, bn, 9.5, MUTED, "middle") + T(xx + (bwid - 12) / 2, cy - ch - 6, f"n={tot:,}", 9, MUTED, "middle")
    out += L(cx, cy, cx + cw, cy, "#8a8578") + T(cx + cw / 2, cy + 32, "family size (distinct loci)", 10.5, MUTED, "middle")
    out += R(cx, cy + 58, 12, 10, GREY) + T(cx + 16, cy + 67, "clustered (one chromosome, < 5 Mb)", 9.5, MUTED)
    out += R(cx + 200, cy + 58, 12, 10, GREY) + R(cx + 200, cy + 58, 12, 10, "url(#hatchC)") + T(cx + 216, cy + 67, "dispersed", 9.5, MUTED)
    out += T(cx, cy + 49, "largest: OR7E 42, GOLGA2/6 37,", 9.5, MUTED) + T(cx + 166, cy + 49, "USP17L 34, SNAR 25", 9.5, ORANGE)
    # D segments vs independent
    dx, dy, dh = 540, 610, 200
    out += panel_label(470, 360, "D", "C1 segments and independent promoters (Σ weights)")
    eff = D["fig2_eff"]
    out += axis_y(dx, dy, dy - dh, 12000, [0, 4000, 8000, 12000], lambda v: f"{v/1000:.0f}k", "promoters")
    for i, (k, col, lab) in enumerate(cls):
        seg, ind = eff[f"C1|{k}"]
        xx = dx + 30 + i * 105
        out += R(xx, dy - seg / 12000 * dh, 34, seg / 12000 * dh, "#ffffff", col, 1.5) + R(xx + 38, dy - ind / 12000 * dh, 34, ind / 12000 * dh, col)
        out += T(xx + 36, dy + 16, lab, 10.5, INK, "middle") + T(xx + 17, dy - seg / 12000 * dh - 5, f"{seg:,}", 9, MUTED, "middle") + T(xx + 55, dy - ind / 12000 * dh - 5, f"{ind:,.0f}", 9, MUTED, "middle")
    out += L(dx, dy, dx + 340, dy, "#8a8578")
    out += R(dx + 230, 380, 10, 10, "#ffffff", INK, 1) + T(dx + 245, 389, "segments", 10, MUTED) + R(dx + 230, 396, 10, 10, INK) + T(dx + 245, 405, "independent", 10, MUTED)
    body = header(2, "Fig 2 · Cohorts and duplicated families",
                  "Promoters are kept in rule-defined cohorts, and duplicated families are detected reliably down to 80% identity with no false links, so 27,833 full-length segments count as 25,669 independent promoters.",
                  "Data ready · to draw")
    body += svg(860, out, "Cohort sizes, family calibration, family sizes and independent promoters")
    body += footer("datasets/cohorts/hpsad_anchors_v7.tsv<br>datasets/families/hpsad_families_v7_summary.json",
                   ["Show C3 and C3x in A as their own bars (alternative anchors) or as a note?", "C: log y (as sketched) keeps the 42-member families visible."])
    return page("Fig 2 sketch", 1240, body)


# ───────────────────────── Fig 3 ─────────────────────────
def fig3():
    out = ""
    prof = D["fig3_prof"]
    specs = [("TATAAA", "TATAAA", 29, BLUE), ("GGGCGG", "GGGCGG / CCGCCC", 56, ORANGE), ("CCAAT", "CCAAT / ATTGG", 74, GREEN)]
    out += panel_label(20, 30, "A", "weighted count per upstream position in C1 (−200 to −1)")
    for i, (k, lab, pos, col) in enumerate(specs):
        x0, y0, pw, ph = 70 + i * 290, 230, 240, 160
        v = prof[k][1:]
        vmax = max(v) * 1.1
        out += T(x0, y0 - ph - 12, lab, 11.5, INK, weight="600") + T(x0 + pw, y0 - ph - 12, f"best k-mer peak −{pos}", 10, MUTED, "end")
        xpos = lambda d: x0 + pw - (d - 1) / 199 * pw
        out += R(xpos(pos + 5), y0 - ph, xpos(pos - 5) - xpos(pos + 5), ph, "#f1eee6")
        pts = " ".join(f"{xpos(d):.1f},{y0 - v[d-1]/vmax*ph:.1f}" for d in range(1, 201))
        out += f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.4"></polyline>'
        base = sum(v[99:199]) / 100
        out += L(x0, y0 - base / vmax * ph, x0 + pw, y0 - base / vmax * ph, MUTED, 1, "4 3")
        out += L(x0, y0, x0 + pw, y0, "#8a8578")
        for d in (200, 100, 1):
            out += T(xpos(d), y0 + 15, f"−{d}", 9.5, MUTED, "middle")
        out += T(x0 - 6, y0 - ph + 8, f"{vmax/1.1:.0f}", 9.5, MUTED, "end") + T(x0 - 6, y0, "0", 9.5, MUTED, "end")
    out += T(70, 268, "shaded = scoring window around the known position · dashed = flank level (mean count over −200 to −101)", 10, MUTED)

    def grid(x0, y0, g, fracs, jits, lab, title):
        o = panel_label(x0 - 50, y0 - 40, lab, title)
        for i, f in enumerate(fracs):
            o += T(x0 - 8, y0 + 22 + i * 30, f"{f*100:g}%", 10, MUTED, "end")
            for j, jt in enumerate(jits):
                n = g.get(f"{f}|{jt}", 0)
                cx = x0 + 30 + j * 70
                o += R(cx - 22, y0 + 8 + i * 30, 44, 22, BLUE if n == 5 else ("#9db5d4" if n else "#f1eee6"), LINE, 0.5, 3)
                o += T(cx, y0 + 23 + i * 30, f"{n}/5", 10.5, "#ffffff" if n == 5 else INK, "middle")
        for j, jt in enumerate(jits):
            o += T(x0 + 30 + j * 70, y0 + 2, f"SD {jt}", 10, MUTED, "middle")
        o += T(x0 - 44, y0 + 2, "carriers", 10, MUTED)
        return o
    out += grid(110, 340, D["fig3_A"], [0.005, 0.01, 0.02, 0.05], [0, 5, 15], "B", "round A: planted motifs recovered")
    out += grid(430, 340, D["fig3_B"], [0.0005, 0.001, 0.002, 0.005], [0, 15, 30], "C", "round B: harder (lower fidelity, fewer carriers)")
    hx, hy = 720, 380
    out += panel_label(690, 300, "D", "round-B k-mers passing FDR")
    tot = D["fig3_hits"]["bench_hits"]
    un = D["fig3_hits"]["bench_unexplained_hits"]
    out += R(hx, hy, 40, 120, BLUE) + R(hx, hy - 2, 40, 2, ORANGE)
    out += T(hx + 50, hy + 60, f"{tot - un:,} planted or real", 10.5, INK) + T(hx + 50, hy + 4, f"{un} neither", 10.5, ORANGE)
    body = header(3, "Fig 3 · A positional statistic that finds known elements and planted motifs",
                  "Scoring each word against its own flanks recovers the TATA box at −29, the GC box at −56 and the CCAAT box at −74 in C1, and every planted motif at fixed or jittered positions.",
                  "Data ready · to draw")
    body += svg(520, out, "Positional profiles and bench recovery")
    body += footer("datasets/discovery/hpsad_bench_summary.json<br>hpsad_bench_hard_summary.json<br>hpsad_discovery_C1_hits.tsv<br>profiles: hpsad_segments_v7.tsv.gz + weights",
                   ["Merge with Fig 4 (motif-level limits) and move this word-level bench to supplementary?", "A: show the expectation the statistic actually uses (own flanks per window) rather than the sketch mean."])
    return page("Fig 3 sketch", 900, body)


# ───────────────────────── Fig 4 ─────────────────────────
def fig4():
    out = ""
    g = D["fig4"]
    carriers = [14, 28, 56, 139, 278, 557, 1392]
    cols = [(0.75, 0), (0.75, 15), (0.85, 0), (0.85, 15), (0.95, 0), (0.95, 15)]
    for p, zone in enumerate(["core", "distal"]):
        x0, y0 = 120 + p * 330, 80
        out += panel_label(x0 - 100, 30, "AB"[p], f"{zone} zone: fraction recovered")
        for j, (f, jt) in enumerate(cols):
            out += T(x0 + 18 + j * 40, y0 - 18, f"{f}", 9.5, MUTED, "middle") + T(x0 + 18 + j * 40, y0 - 6, f"±{jt}", 9.5, MUTED, "middle")
        for i, c in enumerate(carriers):
            out += T(x0 - 8, y0 + 18 + i * 28, f"{c:,}", 10, MUTED, "end")
            for j, (f, jt) in enumerate(cols):
                k = f"{zone}|{c}|{f}|{jt}"
                r, n = g.get(k, [0, 0])
                fr = r / n if n else 0
                shade = f"rgb({int(247 - fr*184)},{int(245 - fr*134)},{int(239 - fr*71)})"
                out += R(x0 + j * 40, y0 + i * 28, 36, 24, shade, LINE, 0.5, 2)
                out += T(x0 + 18 + j * 40, y0 + 16 + i * 28, f"{r}/{n}" if n else "–", 9.5, "#ffffff" if fr > 0.6 else INK, "middle")
        out += T(x0 - 60, y0 - 12, "fidelity", 9.5, MUTED) + T(x0 - 60, y0, "jitter SD", 9.5, MUTED)
        out += f'<text x="{x0 - 72}" y="{y0 + 110}" font-size="10" fill="{MUTED}" text-anchor="middle" transform="rotate(-90 {x0 - 72} {y0 + 110})">carriers</text>'
    bx, by = 760, 300
    out += panel_label(700, 30, "C", "motifs passing FDR, planted cohorts")
    for i, f in enumerate(["0.75", "0.85", "0.95"]):
        v = D["fig4_false"][f]
        xx, acc = bx + i * 46, 0
        for key, col in (("planted", BLUE), ("also_real", GREY), ("other", ORANGE)):
            hh = v[key] / 90 * 200
            out += R(xx, by - acc - hh, 30, hh, col)
            acc += hh
        out += T(xx + 15, by + 14, f, 10, MUTED, "middle") + T(xx + 15, by - acc - 4, str(v["passing"]), 9.5, MUTED, "middle")
    out += L(bx - 6, by, bx + 140, by, "#8a8578") + T(bx + 60, by + 30, "fidelity", 10, MUTED, "middle")
    ly = 345
    for lab, col in (("planted", BLUE), ("also in real data", GREY), ("other", ORANGE)):
        out += R(bx - 6, ly, 10, 10, col) + T(bx + 9, ly + 9, lab, 10, MUTED)
        ly += 15
    body = header(4, "Fig 4 · Detection limits of iterative motif discovery",
                  "A motif carried by about 140 promoters at a fixed position is found; positional jitter of ±15 nt doubles that requirement, and lower fidelity raises it to about 1,400; false hits stay near 3–4%.",
                  "Data ready · to draw")
    body += svg(410, out, "Recovery heatmaps by carriers, fidelity and jitter")
    body += footer("datasets/motifs/hpsad_motifs_bench_zones_f_summary.json<br>(planted into C1; 3 positions per core cell, 2 per distal cell)",
                   ["Colour scale: single-hue blue (as sketched) or viridis?", "Mark the 139-carrier line that the text quotes?"])
    return page("Fig 4 sketch", 790, body)


# ───────────────────────── Fig 5 ─────────────────────────
LETCOL = {"A": GREEN, "C": BLUE, "G": ORANGE, "T": RED}


def logo(x, y, pwm, cw=8, hmax=34):
    o = ""
    for i, col in enumerate(pwm):
        ic = 2 + sum(p * math.log2(p) for p in col if p > 0)
        hs = sorted(((p * ic / 2 * hmax, "ACGT"[k]) for k, p in enumerate(col)), key=lambda t: t[0])
        yb = y
        for h, ch in hs:
            if h < 0.6:
                continue
            sy = h / 7.3
            o += (f'<text x="0" y="0" font-size="10" font-weight="600" fill="{LETCOL[ch]}" text-anchor="middle" '
                  f'transform="translate({x + i*cw + cw/2:.1f},{yb:.1f}) scale({cw/6.2:.2f},{sy:.2f})">{ch}</text>')
            yb -= h
    return o


def fig5_v2():
    ms = D["fig5"]
    groups = [("TATA box (TBP)", ["core:M06", "core:M31", "core:M01"]),
              ("CCAAT box (NF-Y)", ["core:M07", "core:M02", "core:M17"]),
              ("GC box (SP/KLF)", ["core:M10", "core:M04", "core:M19", "core:M27", "core:M13", "core:M15"]),
              ("GC-rich, weak JASPAR match", ["core:M38", "core:M28", "core:M34"]),
              ("CRE (CREB/ATF), incl. half-site at −45", ["core:M12", "core:M09", "core:M39"]),
              ("ETS: class I ×4, class II (CAGGAAGT) ×1", ["core:M11", "core:M03", "core:M05", "core:M33", "core:M26"]),
              ("E-box (TFE3/USF) ×2 · NRF1 · ZBTB33 (Kaiso)", ["core:M08", "core:M14", "core:M20", "core:M25"]),
              ("Unassigned ×4 · last: ZNF-copy artefact", ["core:M24", "core:M23", "core:M18", "core:M40", "core:M16"])]
    by = {m["id"]: m for m in ms}
    out = panel_label(20, 30, "A", "logos of the 32 core motifs passing FDR ≤ 0.05, grouped by element; start position below each")
    y = 52
    for name, ids in groups:
        out += T(30, y + 12, name, 11.5, INK, weight="600")
        x = 30
        for mid in ids:
            m = by[mid]
            out += logo(x, y + 54, m["pwm"], cw=6, hmax=34)
            out += T(x, y + 67, f"−{m['d']}", 10, INK) + T(x, y + 79, m["cons"][:16], 8.5, MUTED, family="IBM Plex Mono, monospace")
            x += max(len(m["pwm"]) * 6, 66) + 22
        y += 90
    out += T(30, y + 4, "letter height = bits (max 2) · pseudogenes and the −17 motif: Fig 7", 10, MUTED)
    # B ratios
    bx, by0 = 70, y + 230
    out += panel_label(20, y + 34, "B", "fold enrichment (O/E) of each motif; window width beside the bar")
    order = sorted(ms, key=lambda m: m["d"])
    vmax = 20
    out += axis_y(bx, by0, by0 - 160, vmax, [0, 5, 10, 15, 20], title="O/E")
    for i, m in enumerate(order):
        xx = bx + 8 + i * 14
        hh = min(m["ratio"], vmax) / vmax * 160
        if m["ratio"] > vmax:
            out += T(xx + 6, by0 - 164, f"{m['ratio']:.0f}↑", 8.5, MUTED, "middle")
        out += R(xx, by0 - hh, 10, hh, BLUE if m["w"] <= 11 else "#9db5d4")
        out += f'<text x="{xx+8:.1f}" y="{by0+8:.1f}" font-size="8.5" fill="{MUTED}" text-anchor="end" transform="rotate(-60 {xx+8:.1f} {by0+8:.1f})">−{m["d"]}</text>'
    out += L(bx, by0, bx + 8 + len(order) * 14, by0, "#8a8578")
    out += R(bx + 300, by0 - 170, 10, 10, BLUE) + T(bx + 315, by0 - 161, "sharp window (≤ 11 nt)", 10, MUTED)
    out += R(bx + 300, by0 - 154, 10, 10, "#9db5d4") + T(bx + 315, by0 - 145, "broad window (≥ 21 nt)", 10, MUTED)
    # C scores
    cx, cy0 = 620, y + 230
    out += panel_label(570, y + 34, "C", "motif scores: real (filled) vs best of each null (open)")
    out += axis_y(cx, cy0, cy0 - 160, 300, [10, 30, 100, 300], title="score (−log10 P)", log=True, vmin=8)
    sc = lambda v: cy0 - (math.log10(max(v, 8)) - math.log10(8)) / (math.log10(300) - math.log10(8)) * 160
    for j, zone in enumerate(["core", "distal"]):
        xx = cx + 50 + j * 110
        for zz, s, p in D["fig5_all"]:
            if zz == zone:
                out += f'<circle cx="{xx - 14 + (int(float(s)*7) % 9):.1f}" cy="{sc(float(s)):.1f}" r="3" fill="{BLUE if p else "#c9c4b6"}"></circle>'
        for s in D["fig5_null"][zone]:
            out += f'<circle cx="{xx + 22:.1f}" cy="{sc(s):.1f}" r="3.5" fill="none" stroke="{ORANGE}" stroke-width="1.5"></circle>'
        out += T(xx + 4, cy0 + 16, zone, 10.5, INK, "middle")
    out += L(cx, cy0, cx + 260, cy0, "#8a8578")
    out += T(cx + 270, cy0 - 150, "blue = passes FDR", 9.5, MUTED, "end") + T(cx + 270, cy0 - 137, "grey = fails · orange = null maxima", 9.5, MUTED, "end")
    h = y + 270
    body = header(5, "Fig 5 · Positioned motifs in full-length promoters of transcribed genes",
                  "In C1 without pseudogenes, all 32 motifs that pass sit within 300 nt of the start site and form the known core-promoter grammar, NRF1 and ZBTB33 sites included, at their known positions; no distal motif beats the shuffled nulls.",
                  "Data ready · to draw")
    body += svg(h, out, "Motif logos, enrichment and scores for C1")
    body += footer("datasets/motifs/hpsad_motifs_C1_transcribed_zones_v2.tsv (pwm_ACGT)<br>hpsad_elements_C1_transcribed_v2.tsv<br>hpsad_motif_annotation_C1_transcribed_v2.tsv<br>…_C1_transcribed_zones_v2_summary.json (null maxima)",
                   ["A: one logo per element (best motif) in the main figure, all 28 in S1 Table?", "Add JASPAR logos beside each element for comparison?"])
    return page("Fig 5 sketch", h + 390, body)


# ───────────────────────── Fig 5 (v4: C1_evidence, three replicates) ─────────────────────────
EL5 = {"WWWAAARN": ("TATA · AT-rich", RED), "YATAAAN": ("TATA box", RED),
       "SNSCCGGAAGTGRCSN": ("ETS class I", GREEN), "NATTTCC": ("ETS-like GGAAATG", "#7aa98a"),
       "SCSGTGACGTSACS": ("CRE", ORANGE), "SSSSTSACGTSACSSS": ("CRE · ACGT palindrome", ORANGE),
       "SSNNGGGGCGGGGC": ("GC box", BLUE), "SGCCCCGCCCCSNS": ("GC box", BLUE),
       "RRCCAATSRGNR": ("CCAAT box", "#7a5ea8"), "NCYGATTGGCYS": ("CCAAT box (reverse)", "#7a5ea8"),
       "NNGRACTACAWYTCC": ("THAP11 / Ronin", "#8a6d3b")}
EL5B = {"NSNSTATAAAWRS": "TATA box, strict", "GSGSCGGARGYGSSG": "GC-rich ETS-like", "CSSSSCYYCCKSSS": "broad GC-rich",
        "SSNSTGCGCASGCG": "NRF1", "SSGGTCACGTGACSS": "E-box CACGTG", "NCATGTGACN": "E-box-like CATGTGAC",
        "CSNSGCCCCGCCCCN": "GC box", "SSCKCYGATTGGCTGS": "CCAAT box, 2nd site"}


def fig5_v4():
    F = D["fig5v4"]
    rows = F["rows"]
    stable = sorted([r for r in rows if r["npass"] == 3], key=lambda r: -r["d"])
    # ── A: syntax map on one position axis
    x0, x1, dlo, dhi = 190, 900, 101, 17          # upstream −101 (left) … −17 (right)
    ppn = (x1 - x0) / (dlo - dhi)
    X = lambda d: x0 + (dlo - d) * ppn             # left edge of upstream position −d
    out = panel_label(20, 30, "A", "the 11 motifs passing in all three replicates, each placed at its start position")
    y = 58
    rh = 44
    for r in stable:
        name, col = EL5[r["cons"]]
        out += T(30, y + 22, name, 11.5, INK, weight="600") + T(30, y + 35, r["cons"], 9, MUTED, family="IBM Plex Mono, monospace")
        out += R(X(r["dmax"]), y + 32, max((r["dmax"] - r["dmin"]) * ppn, 2), 5, col, rx=2, extra=' opacity="0.35"')
        out += logo(X(r["d"]), y + 30, r["pwm"], cw=ppn, hmax=28)
        out += L(X(r["d"]), y + 30, X(r["d"]), y + 39, col, 2)
        y += rh
    for d in range(100, 16, -10):
        out += L(X(d), 50, X(d), y + 2, "#eeeae0") + T(X(d), y + 16, f"−{d}", 10, MUTED, "middle")
    out += L(x0, y + 2, x1, y + 2, "#8a8578") + T(x1, y + 32, "upstream position (nt) · TSS → right", 10, MUTED, "end")
    out += T(30, y + 32, "logo = best replicate · tick = its start · tinted bar = start range over the 3 replicates", 10, MUTED)
    yA = y + 48
    # ── B: replicate matrix
    out += panel_label(20, yA + 10, "B", "each motif in each replicate (score)")
    by = yA + 34
    groups = [("passes 3/3", [r for r in rows if r["npass"] == 3]), ("passes 2/3", [r for r in rows if r["npass"] == 2]),
              ("found, below cut", [r for r in rows if r["npass"] == 0])]
    for j, t in enumerate(["rep 1", "rep 2", "rep 3"]):
        out += T(300 + j * 44 + 18, by, t, 10, MUTED, "middle")
    by += 8
    for gname, gr in groups:
        out += T(30, by + 12, gname, 10, MUTED, style="italic")
        by += 16
        for r in sorted(gr, key=lambda r: r["d"]):
            nm = EL5.get(r["cons"], (EL5B.get(r["cons"], r["cons"]),))[0]
            out += T(30, by + 11, nm, 10.5) + T(200, by + 11, f"−{r['dmin']}" + (f" to −{r['dmax']}" if r["dmax"] != r["dmin"] else ""), 10, MUTED)
            for j, c in enumerate(r["cells"]):
                cx = 300 + j * 44
                if c is None:
                    out += R(cx, by, 36, 15, "none", "#e3dfd3", 1, 3)
                elif c[1]:
                    out += R(cx, by, 36, 15, BLUE, rx=3) + T(cx + 18, by + 11, c[0], 9.5, "#ffffff", "middle")
                else:
                    out += R(cx, by, 36, 15, "#ffffff", BLUE, 1.2, 3) + T(cx + 18, by + 11, c[0], 9.5, BLUE, "middle")
            by += 18
        by += 4
    yB = by
    # ── C: O/E of the stable motifs
    cx0, cy = 520, yA + 34
    out += panel_label(480, yA + 10, "C", "fold enrichment (O/E), mean of 3 replicates")
    vmax = 12
    for r in sorted(stable, key=lambda r: -r["oe"]):
        name, col = EL5[r["cons"]]
        out += T(cx0 + 120, cy + 11, name, 10.5, INK, "end")
        out += R(cx0 + 128, cy + 2, r["oe"] / vmax * 240, 11, col, rx=2)
        out += T(cx0 + 134 + r["oe"] / vmax * 240, cy + 11, f"{r['oe']:.1f}  w{r['w']}", 9.5, MUTED)
        cy += 18
    out += L(cx0 + 128, yA + 30, cx0 + 128, cy, "#8a8578")
    yC = cy
    # ── D: scores vs nulls per replicate
    yD0 = max(yB, yC) + 30
    out += panel_label(20, yD0, "D", "scores: passing (blue), failing (grey), null maxima (orange)")
    dx, dy0, dh = 80, yD0 + 200, 165
    out += axis_y(dx, dy0, dy0 - dh, 300, [10, 30, 100, 300], title="score (−log10 P)", log=True, vmin=8)
    sc = lambda v: dy0 - (math.log10(max(v, 8)) - math.log10(8)) / (math.log10(300) - math.log10(8)) * dh
    k = 0
    for t, lab in zip(["_v4s1", "_v4s2", "_v4s3"], ["rep 1", "rep 2", "rep 3"]):
        S = F["summ"][t]
        for zone in ("core", "distal"):
            xx = dx + 30 + k * 62
            for zz, s, p in S["all"]:
                if zz == zone:
                    out += f'<circle cx="{xx - 10 + (int(s*7) % 9):.1f}" cy="{sc(s):.1f}" r="2.6" fill="{BLUE if p else "#c9c4b6"}"></circle>'
            for s in S["null"][zone]:
                out += f'<circle cx="{xx + 16:.1f}" cy="{sc(s):.1f}" r="3" fill="none" stroke="{ORANGE}" stroke-width="1.3"></circle>'
            out += T(xx + 4, dy0 + 14, zone, 9.5, MUTED, "middle")
            k += 1
        out += T(dx + 30 + (k - 1.5) * 62 + 4, dy0 + 28, lab, 10.5, INK, "middle")
    out += L(dx, dy0, dx + 30 + k * 62, dy0, "#8a8578")
    # ── E: coverage
    ex0, ey = 520, yD0 + 26
    out += panel_label(480, yD0, "E", "promoters carrying ≥ 1 specific motif at its position (rep 1)")
    for nm, real, shuf, n in F["cov"]:
        out += T(ex0 + 110, ey + 12, nm, 10.5, INK, "end") + T(ex0 + 110, ey + 25, f"n = {n:,}", 9, MUTED, "end")
        out += R(ex0 + 118, ey + 2, real / 0.4 * 230, 11, BLUE, rx=2) + T(ex0 + 124 + real / 0.4 * 230, ey + 12, pct(real), 9.5, INK)
        out += R(ex0 + 118, ey + 16, shuf / 0.4 * 230, 9, "#ffffff", "#8a8578", 1, 2) + T(ex0 + 124 + shuf / 0.4 * 230, ey + 24, pct(shuf) + " shuffled", 9, MUTED)
        ey += 34
    out += L(ex0 + 118, yD0 + 22, ex0 + 118, ey, "#8a8578")
    h = max(dy0 + 40, ey + 10)
    body = header(5, "Fig 5 · Eleven reproducible motifs lay out the core-promoter grammar between −29 and −93",
                  "In 12,084 full-length promoters with an EPDnew or CAGE start site, 11 motifs pass in three independent replicates: TATA at −29/−31, ETS, CRE and GC box at −38 to −71, CCAAT at −70 to −82 and THAP11/Ronin at −77 to −93. NRF1 and the E-box pass in two of three; nulls now peak at 26–34.",
                  "Data ready (v4, 13 Oct) · to draw")
    body += svg(h, out, "Syntax map, replicate matrix, enrichment, scores and coverage for C1 with start-site evidence")
    body += footer("datasets/motifs/hpsad_motifs_C1_evidence_zones_v4s1–3.tsv (pwm_ACGT)<br>hpsad_stability_C1_evidence.tsv<br>hpsad_motif_annotation_C1_evidence_v4s1–3.tsv<br>…_v4s1–3_summary.json (null maxima)<br>datasets/coverage/coverage_C1_evidence_v4s1_summary.json",
                   ["A: one row per motif (as sketched) or one row per element, with both orientations of CCAAT and both CRE motifs merged?",
                    "Add the second CCAAT box at −120 to A from a positional profile (planned) rather than from discovery?",
                    "E: coverage from replicate 1, or recompute from the 11 stable motifs only?"])
    return page("Fig 5 sketch", h + 270, body)


# ───────────────────────── Fig 5 (v5: positional profiles of reproducible elements) ─────────────────────────
ORDER5 = ["TATA, AT-rich", "TATA box", "ETS class I", "CRE (CREB/ATF)", "GC box (SP/KLF)", "GCGC repeat (ZNF131)",
          "GGAGGAG (ZNF263/MAZ)", "AP-1-like TGAGTCA", "CCAAT box (NF-Y)", "THAP11 / Ronin"]


def fig5():
    F = D["fig5v5"]
    rows = {r["name"]: r for r in F["rows"]}
    rows = [rows[n] for n in ORDER5]
    # ── A: logo + profile per element on one position axis
    px0, px1, dlo, dhi = 300, 900, -200, -1
    X = lambda p: px0 + (p - dlo) / (dhi - dlo) * (px1 - px0)
    out = panel_label(20, 30, "A", "each element's enrichment along the promoter (observed / expected; expectation = 10 shuffled cohorts)")
    y0, rh, ph = 52, 50, 36
    for r in rows:
        out += T(30, y0 + 15, r["name"], 11, INK, weight="600")
        out += T(30, y0 + 28, r["cons"], 8.5, MUTED, family="IBM Plex Mono, monospace")
        out += logo(180, y0 + 42, r["pwm"], cw=6.5, hmax=30)
        ymax = max(3.0, max(r["y"]) * 1.05)
        base = y0 + 44
        Y = lambda v: base - (min(v, ymax) - 1) / (ymax - 1) * ph if v >= 1 else base
        pts = " ".join(f"{X(p):.1f},{Y(v):.1f}" for p, v in zip(r["pos"], r["y"]))
        out += (f'<polygon points="{X(r["pos"][0]):.1f},{base} {pts} {X(r["pos"][-1]):.1f},{base}" fill="{BLUE}" '
                f'fill-opacity="0.85" stroke="#2c5282" stroke-width="0.6"></polygon>')
        out += L(px0, base, px1, base, "#bdb8aa", 0.8)
        out += T(px0 - 6, base - ph + 4, f"{ymax:.0f}×", 8.5, MUTED, "end") + T(px0 - 6, base + 3, "1", 8.5, MUTED, "end")
        for k, pk in enumerate(r["peaks"][:3]):
            if dlo <= pk["c"] <= dhi:
                out += (f'<polygon points="{X(pk["c"]) - 3.5:.1f},{base + 9} {X(pk["c"]) + 3.5:.1f},{base + 9} {X(pk["c"]):.1f},{base + 3}" '
                        f'fill="{ORANGE if k == 0 else "#e8b48a"}"></polygon>')
        y0 += rh
    for p in (-200, -175, -150, -125, -100, -75, -50, -25):
        out += L(X(p), 48, X(p), y0 + 2, "#eeeae0", 0.8) + T(X(p), y0 + 16, f"−{-p}", 10, MUTED, "middle")
    out += L(X(-40), 48, X(-40), y0 + 2, RED, 1, "4 3")
    out += T(X(-40) + 4, 46, "−40", 9.5, RED)
    out += T(px1, y0 + 32, "motif start, upstream position (nt) · TSS → right", 10, MUTED, "end")
    out += T(30, y0 + 32, "▲ summits (orange = highest) · logo = representative PWM", 10, MUTED)
    yA = y0 + 52
    # ── B: element x replicate
    out += panel_label(20, yA, "B", "result of each element in each replicate (best motif score)")
    by = yA + 24
    cx0 = 210
    out += T(cx0 + 6 * 38 + 62, by - 11, "sense share", 10, MUTED) + T(cx0 + 6 * 38 + 62, by + 1, "obs / exp", 8.5, MUTED)
    out += T(cx0 + 3 * 38 + 9, by, "replicates 1–6", 10, MUTED, "middle")
    by += 8
    for r in rows:
        out += T(30, by + 11, r["name"], 10.5)
        for j, c in enumerate(r["cells"]):
            x = cx0 + j * 38 + 9
            if c is None:
                out += R(x, by, 32, 15, "none", "#e3dfd3", 1, 3)
            elif c[1]:
                out += R(x, by, 32, 15, BLUE, rx=3) + T(x + 16, by + 11, c[0], 9, "#ffffff", "middle")
            else:
                out += R(x, by, 32, 15, "#ffffff", BLUE, 1.2, 3) + T(x + 16, by + 11, c[0], 9, BLUE, "middle")
        n = sum(1 for c in r["cells"] if c and c[1])
        out += T(cx0 + 6 * 38 + 26, by + 11, f"{n}/6", 10, INK if n == 6 else MUTED, weight="600" if n == 6 else None)
        sd = F["strand"][r["cons"]]
        if sd["pal"]:
            txt, col, wt = "palindrome", MUTED, None
        else:
            sig = sd["p"] is not None and sd["p"] < 0.001
            txt = f"{sd['obs']:.2f} / {sd['exp']:.2f}" if sd["exp"] is not None else f"{sd['obs']:.2f} / 0.50"
            col, wt = (ORANGE, "600") if sig else (INK, None)
        out += T(cx0 + 6 * 38 + 62, by + 11, txt, 9.5, col, weight=wt)
        by += 18
    out += T(30, by + 14, "filled = passes FDR · open = found below the cut · blank = not found", 9.5, MUTED)
    out += T(30, by + 28, "sense share at the main summit; orange = differs from expected (binomial p < 0.001)", 9.5, MUTED)
    by += 14
    # ── C: summit scores vs null summits
    cx, cw_, cy = 745, 160, yA + 60
    out += panel_label(590, yA, "C", "top summit score vs null profiles")
    sx = lambda v: cx + (math.log10(max(v, 1)) / math.log10(300)) * cw_
    for v in F["nulls"]:
        out += L(sx(v), cy - 6, sx(v), cy + 6, "#b7b2a4", 1)
    out += T(cx - 10, cy + 4, f"null profiles (n = {len(F['nulls'])})", 9.5, MUTED, "end")
    out += L(sx(F["cut"]), cy - 18, sx(F["cut"]), cy + 26 + len(rows) * 15, RED, 1, "4 3")
    out += T(sx(F["cut"]) + 4, cy + 22, "cut", 9.5, RED)
    yy = cy + 34
    for r in sorted(rows, key=lambda r: -max(p["s"] for p in r["peaks"])):
        s = max(p["s"] for p in r["peaks"])
        out += f'<circle cx="{sx(s):.1f}" cy="{yy:.1f}" r="3.5" fill="{BLUE}"></circle>'
        out += T(cx - 10, yy + 3.5, r["name"], 9.5, INK, "end") + L(cx - 6, yy, sx(s) - 5, yy, "#eeeae0", 0.8)
        yy += 15
    out += L(cx, yy + 4, cx + cw_, yy + 4, "#8a8578")
    for v in (1, 3, 10, 30, 100, 300):
        out += L(sx(v), yy + 4, sx(v), yy + 8, "#8a8578") + T(sx(v), yy + 19, v, 9.5, MUTED, "middle")
    out += T(cx + cw_, yy + 33, "summit score (−log10 P), log scale", 9.5, MUTED, "end")
    # ── D: TATA+ vs TATA- profiles
    yD = max(by + 30, yy + 44) + 20
    nT = int(next(r for r in F["slots"] if r["group"] == "TATA+")["promoters_in_group"])
    nN = int(next(r for r in F["slots"] if r["group"] == "TATA-")["promoters_in_group"])
    out += panel_label(20, yD, "D", f"profiles in TATA promoters (orange, n = {nT:,}) and TATA-less promoters (blue, n = {nN:,})")
    dx0, dx1, dl, dh_ = 300, 900, -150, -1
    XD = lambda p: dx0 + (p - dl) / (dh_ - dl) * (dx1 - dx0)
    rowsD = [("YATAAA", "TATA box"), ("STGACGTCAS", "CRE (CREB/ATF)"), ("SNSGGGGCGGGGCS", "GC box (SP/KLF)"),
             ("SNNSCCGGAAGTGRCS", "ETS class I"), ("RGCCAATSRG", "CCAAT box (NF-Y)")]
    y = yD + 22
    out += R(XD(-58), y - 4, XD(-46) - XD(-58), len(rowsD) * 70, "#eeeae0") + R(XD(-37), y - 4, XD(-18) - XD(-37), len(rowsD) * 70, "#f6e3d4")
    out += T((XD(-58) + XD(-46)) / 2, y - 8, "slot A", 9, MUTED, "middle") + T((XD(-37) + XD(-18)) / 2, y - 8, "TATA position", 9, MUTED, "middle")
    for cons, nm in rowsD:
        ymx = max(3.0, max(max(F["tata_prof"][f"{cons}|{g}"]["y"]) for g in ("TATA+", "TATA-")) * 1.05)
        base = y + 60
        YD = lambda v: base - (min(max(v, 0), ymx)) / ymx * 54
        out += T(30, y + 34, nm, 11, INK, weight="600")
        for g, col in (("TATA-", BLUE), ("TATA+", ORANGE)):
            pr = F["tata_prof"][f"{cons}|{g}"]
            pts = " ".join(f"{XD(p):.1f},{YD(v):.1f}" for p, v in zip(pr["pos"], pr["y"]))
            out += f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.4"></polyline>'
        out += L(dx0, YD(1), dx1, YD(1), "#bdb8aa", 0.6, "2 2") + L(dx0, base, dx1, base, "#8a8578", 0.6)
        out += T(dx0 - 6, base - 50, f"{ymx:.0f}×", 8.5, MUTED, "end") + T(dx0 - 6, YD(1) + 3, "1", 8.5, MUTED, "end")
        y += 70
    for p in (-150, -125, -100, -75, -50, -25):
        out += T(XD(p), y + 10, f"−{-p}", 10, MUTED, "middle")
    out += T(dx1, y + 26, "motif start, upstream position (nt)", 10, MUTED, "end")
    h = y + 40
    body = header(5, "Fig 5 · Reproducible promoter elements and where they sit",
                  "Seven elements pass in all six replicates. TATA is pinned at −31 and oriented; CRE, GC box, ETS, GGAGGAG and AP-1-like summits sit at −53 to −55 in all promoters; in TATA-less promoters ETS, CRE and CCAAT also take the TATA position (−18 to −37) — one slot, two kinds of occupant.",
                  "Final weights (_f1–_f6, 17 Oct) · to draw")
    body += svg(h, out, "Positional profiles, replicate table and null calibration for reproducible promoter elements")
    body += footer("datasets/profiles/profile_C1_evidence_f1_6rep.tsv (O, E per position)<br>syntax_strand_ · syntax_tata_profile_ · syntax_tata_slot_C1_evidence_f1_6rep.tsv<br>profile_peaks_C1_evidence_f1_6rep.tsv (summits)<br>profile_C1_evidence_f1_6rep_summary.json (null summits, cut)<br>datasets/motifs/hpsad_motifs_C1_evidence_zones_f1–f6.tsv",
                   ["D: keep all six rows, or only CRE, GC box and ETS for the main figure?",
                    "Supplementary Fig S2: the motif-level replicate matrix (previous version of this board)",
                    "Coverage moves to the promoter-families figure"])
    return page("Fig 5 sketch", h + 260, body)


# ───────────────────────── Fig 6 (start-site region) ─────────────────────────
def fig_tss():
    F = D["fig6tss"]
    rows = {r["name"]: r for r in F["rows"]}
    px0, px1, dlo, dhi = 250, 900, -60, 100
    X = lambda p: px0 + (p - dlo) / (dhi - dlo) * (px1 - px0)

    def track(y0, r, col):
        ymax = max(3.0, max(r["y"]) * 1.05)
        base = y0 + 40
        Y = lambda v: base - (min(v, ymax) - 1) / (ymax - 1) * 32 if v >= 1 else base
        pts = " ".join(f"{X(p):.1f},{Y(v):.1f}" for p, v in zip(r["pos"], r["y"]) if dlo <= p <= dhi)
        ps = [p for p in r["pos"] if dlo <= p <= dhi]
        o = (f'<polygon points="{X(ps[0]):.1f},{base} {pts} {X(ps[-1]):.1f},{base}" fill="{col}" fill-opacity="0.85" '
             f'stroke="{col}" stroke-width="0.6"></polygon>')
        o += L(px0, base, px1, base, "#bdb8aa", 0.8)
        o += T(px0 - 6, base - 28, f"{ymax:.0f}×", 8.5, MUTED, "end") + T(px0 - 6, base + 3, "1", 8.5, MUTED, "end")
        o += T(30, y0 + 20, r["name"], 11, INK, weight="600") + T(30, y0 + 33, r["cons"], 8.5, MUTED, family="IBM Plex Mono, monospace")
        return o

    out = panel_label(20, 30, "A", "core-promoter elements around the start site (observed / expected, sense strand)")
    y = 50
    order = ["TATA box", "Initiator, CA+1G", "Initiator, TCA+1TT", "TCT (consensus)", "DPE (consensus)", "YY1"]
    for n in order:
        out += track(y, rows[n], BLUE)
        y += 50
    yA_end = y
    out += panel_label(20, y + 16, "B", "transcript features in the same window (not promoter elements)")
    y += 30
    for n in ["start codon (Kozak)", "5' splice donor"]:
        out += track(y, rows[n], "#9a958a")
        y += 50
    for p in (-60, -40, -20, 20, 40, 60, 80, 100):
        out += L(X(p), 46, X(p), y + 2, "#eeeae0", 0.8) + T(X(p), y + 16, f"{p:+d}", 10, MUTED, "middle")
    out += L(X(0), 46, X(0), yA_end + 4, RED, 1, "4 3") + L(X(0), yA_end + 38, X(0), y + 2, RED, 1, "4 3") + T(X(0) + 4, 44, "start site", 9.5, RED)
    for p, lab in ((28, "DPE +28"), (-31, "TATA −31")):
        out += T(X(p), 44, lab, 8.5, MUTED, "middle")
    out += T(px1, y + 32, "motif start relative to the start site (nt; +1 = start-site base)", 10, MUTED, "end")
    # ── C: targeted test
    yC = y + 56
    out += panel_label(20, yC, "C", "published consensus at its published position: O/E (promoters matching)")
    groups = ["all", "TATA+", "TATA- CpG", "TATA- non-CpG", "Initiator+ TATA-"]
    gl = ["all", "TATA+", "TATA− CpG", "TATA− non-CpG", "Inr+ TATA−"]
    names = []
    for r in F["targeted"]:
        if r["element"] not in names:
            names.append(r["element"])
    cx0, cw = 330, 112
    for j, g in enumerate(gl):
        out += T(cx0 + j * cw + cw / 2, yC + 24, g, 10, MUTED, "middle")
    yy = yC + 32
    for nm in names:
        rr = {r["group"]: r for r in F["targeted"] if r["element"] == nm}
        r0 = next(iter(rr.values()))
        out += T(30, yy + 12, nm.replace("(control)", "").strip(), 10.5, INK) + T(250, yy + 12, f"{int(r0['start_from']):+d}..{int(r0['start_to']):+d}", 9, MUTED)
        for j, g in enumerate(groups):
            r = rr[g]
            ratio, sc, n = float(r["ratio"]), float(r["score"]), int(r["promoters_matching"])
            few = n < 5 or float(r["E"]) < 1
            sig = sc >= 3 and ratio > 1 and not few
            fill = BLUE if sig and ratio >= 3 else ("#7d9fc7" if sig and ratio >= 1.5 else ("#c3d3e6" if sig else "#f1eee6"))
            ink = "#ffffff" if sig and ratio >= 1.5 else INK
            txt = f"{ratio:.0f}" if ratio >= 10 else f"{ratio:.2f}"
            out += R(cx0 + j * cw + 4, yy, cw - 8, 17, fill, rx=3)
            out += T(cx0 + j * cw + cw / 2, yy + 12, "—" if n == 0 else (f"n = {n}" if few else f"{txt}  ({n:,})"), 9.5, ink, "middle")
        if "control" in nm:
            out += T(cx0 + 5 * cw + 6, yy + 12, "control", 9, MUTED)
        yy += 20
    out += T(30, yy + 14, "shaded = more than expected (Poisson −log10 P ≥ 3): dark ≥ 3-fold, mid ≥ 1.5-fold, light < 1.5-fold; beige = not enriched; n = fewer than 5 matching promoters", 9.5, MUTED)
    # ── D: TCT carriers
    yD = yy + 44
    out += panel_label(20, yD, "D", "promoters with a TCT at the start site (64): gene classes")
    groups_txt = [("ribosomal proteins", "RPL5, RPL13, RPL22, RPL30, RPL38, RPLP1, RPS3, RPS16, MRPS25"),
                  ("translation factors", "EEF1A1, EIF2D, EIF3A, EIF3E, NACA2"),
                  ("snoRNA hosts / RNA processing", "SNHG12, SNHG16, SNHG26, FBL, NONO"),
                  ("other", "COX7C, CDK7, SEC11A, TOP1MT and 41 more (11 lncRNAs or uncharacterised loci)")]
    yt = yD + 22
    for a_, b_ in groups_txt:
        out += T(30, yt, a_, 10.5, INK, weight="600") + T(250, yt, b_, 10.5, INK)
        yt += 17
    h = yt + 16
    body = header(6, "Fig 6 · The start-site region: Initiator, TCT, DPE and YY1 at fixed positions",
                  "Extending each promoter to +100 adds the Initiator (CA at −1/+1, 2,730 promoters), TCT at the start site in ribosomal-protein and translation genes, a weak DPE at +28 only in TATA-less promoters, and YY1 at about +10 and +21. MTE, BREu, BREd and XCPE1 show no positional signal.",
                  "Final weights (_f, 17 Oct) · to draw")
    body += svg(h, out, "Start-site region profiles, targeted consensus test and TCT carriers")
    body += footer("datasets/profiles/downstream_profiles_C1ev_down_f.tsv<br>downstream_targeted_C1ev_down_f.tsv<br>datasets/motifs/hpsad_motifs_C1ev_down*_zones_f.tsv<br>composite sequences: hpsad_segments_v7_down100.tsv.gz",
                   ["A: DPE consensus also peaks at +1 because RGWYVT overlaps the Initiator's AG; keep, or mask the Initiator for this row?",
                    "B: keep the transcript features here or move them to a supplementary figure?"])
    return page("Fig 6 sketch", h + 270, body)


# ───────────────────────── Fig 6 ─────────────────────────
def fig6():
    F = json.load(open("fig7_f.json"))
    runs = ["C1 Pol II", "C1 coding", "C1 lncRNA", "C1 small RNA", "C1 pseudo", "C2", "C2b", "C3", "C4 Pol II", "C4 small RNA"]
    short = ["C1 Pol II", "coding", "lncRNA", "small RNA", "pseudo", "C2", "C2b", "C3", "C4 Pol II", "C4 sm. RNA"]
    els = [("TATA box", BLUE), ("Initiator", BLUE), ("YY1", BLUE), ("CCAAT box", BLUE), ("GC box", BLUE), ("ETS", BLUE),
           ("CRE", BLUE), ("NRF1", BLUE), ("E-box", BLUE), ("THAP11/Ronin", BLUE),
           ("lncRNA CRE-like", "#7a62b3"), ("CCTTTAAGG (ZNF536-like)", "#7a62b3"),
           ("Start codon context", GREY), ("5' splice donor", GREY),
           ("tRNA B box", GREEN), ("C/D box C (snoRNA)", GREEN)]
    out = panel_label(20, 28, "A", "elements passing FDR ≤ 0.05 on −300..+100, by run (box: start of best motif; area = carriers, shade = O/E)")
    x0, y0, cw, rh = 230, 78, 68, 27
    for j, r in enumerate(short):
        out += T(x0 + j * cw + cw / 2, y0 - 12, r, 10, INK, "middle", weight="600" if j == 0 else None)
    for i, (nm, col) in enumerate(els):
        yy = y0 + i * rh
        if i % 2 == 0:
            out += R(20, yy - 4, x0 + 10 * cw - 20, rh, "#f4f1e8")
        out += T(30, yy + 13, nm.replace("'", "′"), 11, col if col != GREY else MUTED)
        for j, run in enumerate(runs):
            e = F[run]["elements"].get(nm)
            cx = x0 + j * cw + cw / 2
            if not e:
                out += f'<circle cx="{cx:.1f}" cy="{yy + 9:.1f}" r="1.6" fill="#b9b4a6"></circle>'
                continue
            b = max(e, key=lambda m: m["score"])
            side = 8 + min(1.0, math.sqrt(b["n"] / 2000)) * 16
            op = 0.25 + 0.75 * min(1.0, math.log10(max(b["ratio"], 1.01)) / 2)
            out += R(cx - side / 2, yy + 9 - side / 2, side, side, col, rx=2, extra=f' fill-opacity="{op:.2f}"')
            out += T(cx + side / 2 + 2, yy + 13, f"{b['pos']:+d}", 8.5, INK)
            if len(e) > 1:
                out += T(cx - side / 2 - 2, yy + 13, str(len(e)), 8, MUTED, "end")
    yb = y0 + len(els) * rh + 6
    out += T(x0 - 8, yb + 10, "passing motifs", 10, MUTED, "end")
    out += T(x0 - 8, yb + 26, "unassigned, O/E>20, <150 carriers", 10, MUTED, "end")
    for j, run in enumerate(runs):
        cx = x0 + j * cw + cw / 2
        u = sum(1 for m in F[run]["unassigned"] if m["ratio"] > 20 and m["n"] < 150)
        out += T(cx, yb + 10, str(F[run]["n_pass"]), 10, INK, "middle") + T(cx, yb + 26, str(u), 10, ORANGE if u > 5 else INK, "middle")
    net, bottom = arch_net(yb + 62)
    out += net
    body = header(7, "Fig 7 · The core grammar recurs in every Pol II cohort and combines into a small set of architectures with a TATA pole and a CpG pole",
                  "Every Pol II cohort recovers the core elements at the same positions, and small-RNA genes give their own internal elements (A). In the main cohort, 54% of promoters carry an element at its position (9% of shuffled copies), 22 architectures hold 97% of promoters, and they sort into an AT-rich TATA pole and a CpG-island pole (B).",
                  "Final weights (_f, 17 Oct) · to draw")
    body += svg(bottom, out, "Element presence by cohort and promoter architectures")
    body += footer("A: datasets/motifs/hpsad_motifs_RUN_zones_f.tsv, hpsad_elements_RUN_f.tsv; fig7_f.json<br>B: datasets/profiles/architecture_nodes_C1ev_v7.tsv, architecture_summary_C1ev_v7.json<br>families: hpsad_families_v7_down100.tsv",
                   ["Decided 16 Oct: old B (few-carrier scatter) folded into the counts row of A; architecture network is the new B; element call rates go to S1 Fig.", "Weak calls (GC box) kept with dashed outlines; no edges from 'no element'."])
    return page("Fig 7 sketch", bottom + 260, body)


# ───────────────────────── Fig 7 ─────────────────────────
def fig7():
    out = panel_label(20, 28, "A", "a copy inserted by LINE-1 (schematic, not to scale)")
    y = 80
    out += L(40, y, 880, y, "#8a8578", 2)
    out += R(150, y - 12, 90, 24, ORANGE, rx=3, extra=' fill-opacity="0.85"') + T(195, y + 4, "TSD", 11, "#ffffff", "middle", "600")
    out += T(150, y - 20, "TT|AAAA…", 11, INK, family="IBM Plex Mono, monospace") + T(195, y + 32, "starts ≈ −16, 12–15 nt", 10, MUTED, "middle")
    out += f'<polygon points="250,{y-14} 600,{y-14} 625,{y} 600,{y+14} 250,{y+14}" fill="{BLUE}" fill-opacity="0.85"></polygon>'
    out += T(425, y + 4, "gene copy (processed pseudogene)", 11, "#ffffff", "middle", "600")
    out += L(250, y - 18, 250, y - 30, MUTED, 1) + T(254, y - 22, "5′ edge", 10, MUTED)
    out += T(663, y - 6, "AAAAAA", 11, BLUE, "middle", family="IBM Plex Mono, monospace") + T(663, y + 32, "poly(A)", 10, MUTED, "middle")
    out += R(700, y - 12, 90, 24, ORANGE, rx=3, extra=' fill-opacity="0.85"') + T(745, y + 4, "TSD", 11, "#ffffff", "middle", "600")
    out += f'<path d="M195,{y+40} C195,{y+70} 745,{y+70} 745,{y+40}" fill="none" stroke="{ORANGE}" stroke-width="1.2" stroke-dasharray="4 3"></path>'
    out += T(470, y + 78, "same sequence", 10, ORANGE, "middle")
    l1 = {(r["cohort"], r["group"]): r for r in D["fig7_l1"]}
    ts = {(r["cohort"], r["group"]): r for r in D["fig7_tsd"]}
    groups = [("C1", "protein_coding", "C1 coding"), ("C1", "non_coding", "C1 non-cod."), ("C1", "pseudogene", "C1 pseudo"),
              ("C1", "pseudogene_RPL_RPS", "C1 RPL/RPS"), ("C4", "pseudogene", "C4 pseudo"), ("C4", "pseudogene_RPL_RPS", "C4 RPL/RPS")]
    # B
    bx, by, bh = 70, 380, 150
    out += panel_label(20, 200, "B", "TTAAAR/TTRAAA at −25 to −10 (filled), background (open), TATAWA at −35 to −25 (hatched)")
    out += '<defs><pattern id="hatch" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="5" stroke="#3f6fa8" stroke-width="2"></line></pattern></defs>'
    out += axis_y(bx, by, by - bh, 0.2, [0, 0.05, 0.1, 0.15, 0.2], lambda v: f"{v*100:.0f}%", "promoters")
    for i, (c, g, lab) in enumerate(groups):
        r = l1[(c, g)]
        xx = bx + 14 + i * 70
        for j, (v, fill, st) in enumerate([(float(r["l1_frac_-25_-10"]), ORANGE, "none"), (float(r["l1_frac_bg"]), "#ffffff", ORANGE), (float(r["tata_frac_-35_-25"]), "url(#hatch)", BLUE)]):
            hh = v / 0.2 * bh
            out += R(xx + j * 17, by - hh, 15, hh, fill, st, 1.2)
        out += T(xx + 25, by + 15, lab, 9.5, INK, "middle")
    out += L(bx, by, bx + 440, by, "#8a8578")
    # C
    cx, cy, ch = 590, 380, 150
    out += panel_label(540, 225, "C", "target-site duplication: real (filled) vs mismatched pairs (open)")
    out += axis_y(cx, cy, cy - ch, 0.16, [0, 0.05, 0.1, 0.15], lambda v: f"{v*100:.0f}%", None)
    for i, (c, g, lab) in enumerate(groups):
        r = ts[(c, g)]
        xx = cx + 10 + i * 50
        v, nl = float(r["tsd_frac"]), float(r["tsd_null"])
        out += R(xx, cy - v / 0.16 * ch, 16, v / 0.16 * ch, ORANGE) + R(xx + 18, cy - max(nl / 0.16 * ch, 1), 16, max(nl / 0.16 * ch, 1), "#ffffff", ORANGE, 1.2)
        out += f'<text x="{xx+18:.1f}" y="{cy+10:.1f}" font-size="9" fill="{INK}" text-anchor="end" transform="rotate(-35 {xx+18:.1f} {cy+10:.1f})">{lab}</text>'
        out += T(xx + 8, cy - v / 0.16 * ch - 4, pct(v), 8.5, MUTED, "middle")
    out += L(cx, cy, cx + 310, cy, "#8a8578")
    # D
    dx, dy, dh = 110, 640, 140
    out += panel_label(20, 440, "D", "pseudogenes: duplication with (filled) and without (open) the −17 word")
    out += axis_y(dx, dy, dy - dh, 0.2, [0, 0.1, 0.2], lambda v: f"{v*100:.0f}%", "with TSD")
    for i, (c, g, lab) in enumerate([("C1", "pseudogene", "C1"), ("C2", "pseudogene", "C2"), ("C4", "pseudogene", "C4"), ("C1", "pseudogene_RPL_RPS", "C1 RPL/RPS")]):
        r = ts[(c, g)]
        xx = dx + 20 + i * 85
        a, b = float(r["tsd_frac_if_scar"]), float(r["tsd_frac_if_no_scar"])
        a = min(a, 0.2)
        out += R(xx, dy - a / 0.2 * dh, 26, a / 0.2 * dh, ORANGE) + R(xx + 30, dy - b / 0.2 * dh, 26, b / 0.2 * dh, "#ffffff", ORANGE, 1.2)
        out += T(xx + 28, dy + 15, lab, 10, INK, "middle")
        if float(r["tsd_frac_if_scar"]) > 0.2:
            out += T(xx + 13, dy - dh - 4, pct(float(r["tsd_frac_if_scar"])) + " ↑", 9, MUTED, "middle")
    out += L(dx, dy, dx + 360, dy, "#8a8578")
    # E
    ex, ey, ew, eh = 560, 640, 320, 140
    out += panel_label(510, 440, "E", "start of the upstream copy (nt before the 5′ edge), pseudogenes")
    dep = {int(k): v for k, v in D["fig7_depth"].items()}
    vmax = 150
    out += R(ex + (15 - 12) / 19 * ew, ey - eh, 5 / 19 * ew, eh, "#f6e7d3")
    out += axis_y(ex, ey, ey - eh, vmax, [0, 50, 100, 150], title=None)
    for d, n in dep.items():
        xx = ex + (d - 12) / 19 * ew
        out += R(xx + 1, ey - n / vmax * eh, ew / 19 - 3, n / vmax * eh, ORANGE)
    for d in (12, 16, 20, 25, 30):
        out += T(ex + (d - 12) / 19 * ew + 7, ey + 15, f"−{d}", 9.5, MUTED, "middle")
    out += L(ex, ey, ex + ew, ey, "#8a8578") + T(ex + ew, ey - eh + 4, "band = −17 motif start", 9.5, MUTED, "end")
    out += T(ex + ew, ey + 32, "−30 bar = search limit (copies that start further up)", 9.5, MUTED, "end")
    body = header(8, "Fig 8 · The −17 motif marks where retrotransposed pseudogenes landed",
                  "The AT-rich word at −17 is enriched only in pseudogenes, avoids TATA boxes, and starts the target-site duplication that flanks retrotransposed copies 60 times more often than chance.",
                  "Data ready (new today) · to draw")
    body += svg(700, out, "LINE-1 schematic, word enrichment and target-site duplications")
    body += footer("datasets/qc/l1_footprint_v7.tsv<br>datasets/qc/tsd_v7_summary.tsv<br>datasets/qc/tsd_v7_genes.tsv (per gene)",
                   ["Keep B in the main figure or move it to supplementary now that C–E show the duplications directly?", "E: add the same histogram for non-coding genes (100 duplications) as a contrast?"])
    return page("Fig 7 sketch", 1090, body)


# ───────────────────────── architecture network (Fig 7B) and element calls (S1 Fig) ─────────────────────────
def arch_colour(f):                      # 0 -> orange (AT-rich, non-CpG), 1 -> blue (CpG island)
    a = (217, 115, 31); b = (63, 111, 168)
    t = max(0.0, min(1.0, (f - 0.1) / 0.7))
    return "#%02x%02x%02x" % tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(3))


def arch_net(y0, letter="B"):
    """architecture tree starting at y0; returns (svg, bottom y)."""
    nodes = list(csv.DictReader(open("architecture_nodes_C1ev_v7.tsv"), delimiter="\t"))
    Sm = json.load(open("architecture_summary_C1ev_v7.json"))
    for n in nodes:
        n["set"] = frozenset() if n["arch"] == "(none)" else frozenset(n["arch"].split(" + "))
        n["w"] = float(n["weighted"]); n["r"] = 5 + math.sqrt(n["w"]) * 0.62
        n["cpg"] = float(n["cpg_fraction"]); n["ratio"] = float(n["ratio"])
    order = ["TATA", "Inr", "TCT", "CCAAT", "GC box", "YY1", "CRE", "ETS", "NRF1", "THAP11", "AP-1"]
    out = panel_label(20, y0, letter, "architectures of C1 promoters with an EPDnew or CAGE start site (node area = promoters; colour = share CpG-island-like)")
    x0, x1 = 70, 880
    by = {k: [n for n in nodes if int(n["n_elements"]) == k] for k in range(4)}
    r1 = max(n["r"] for n in by[1]); r2 = max(n["r"] for n in by[2])
    ys = {1: y0 + 70 + r1, 2: 0, 3: 0}
    ys[2] = ys[1] + r1 + 70 + r2
    band2 = ys[2] + r2 + 13 + 26 + 18                 # below every (staggered) layer-2 label
    ys[3] = band2 + 60
    singles = sorted(by[1], key=lambda n: order.index(next(iter(n["set"]))))
    for i, n in enumerate(singles):
        n["x"] = x0 + 20 + i * (x1 - x0 - 40) / max(1, len(singles) - 1)
    xs = {next(iter(n["set"])): n["x"] for n in singles}
    for k in (2, 3):
        lay = by[k]
        for n in lay:
            n["x"] = sum(xs.get(e, (x0 + x1) / 2) for e in n["set"]) / len(n["set"])
        lay.sort(key=lambda n: n["x"])
        for i in range(1, len(lay)):
            gap = lay[i - 1]["r"] + lay[i]["r"] + 52
            if lay[i]["x"] < lay[i - 1]["x"] + gap:
                lay[i]["x"] = lay[i - 1]["x"] + gap
        shift = max(0, lay[-1]["x"] - x1) if lay else 0
        for n in lay:
            n["x"] -= shift
    for k in (1, 2, 3):
        for n in by[k]:
            n["y"] = ys[k]
    # edges (no edges from 'no element'): layer 1 -> 2 from node bottom to node top (single labels sit above their nodes);
    # layer 2 -> 3 start below the layer-2 label band
    for a_ in by[1] + by[2]:
        for b_ in by[2] + by[3]:
            if len(b_["set"]) == len(a_["set"]) + 1 and a_["set"] < b_["set"]:
                ya = a_["y"] + a_["r"] if len(a_["set"]) == 1 else band2
                out += L(a_["x"], ya, b_["x"], b_["y"] - b_["r"], "#8a8578", 1.2)
                if len(a_["set"]) == 2:
                    out += L(a_["x"], band2 - 6, a_["x"], band2, "#8a8578", 1.2)
    for k in (1, 2, 3):
        for i, n in enumerate(by[k]):
            dash = ' stroke-dasharray="3 2"' if n["ratio"] < 2 else ""
            out += (f'<circle cx="{n["x"]:.1f}" cy="{n["y"]:.1f}" r="{n["r"]:.1f}" fill="{arch_colour(n["cpg"])}" fill-opacity="0.85" '
                    f'stroke="#5b574d" stroke-width="1"{dash}></circle>')
            lab = " + ".join(sorted(n["set"], key=order.index))
            sub = f"{n['w']:,.0f}  (×{n['ratio']:.0f})"
            if k == 1:
                out += T(n["x"], n["y"] - n["r"] - 18, lab, 10, INK, "middle", weight="600") + T(n["x"], n["y"] - n["r"] - 6, sub, 9, MUTED, "middle")
            else:
                ly = n["y"] + n["r"] + 13 + (26 if k == 2 and i % 2 == 1 else 0)
                out += T(n["x"], ly, lab, 10, INK, "middle") + T(n["x"], ly + 12, sub, 9, MUTED, "middle")
    # 'no element' and legend, top right of the panel
    lx, ly0 = 690, y0 + 2
    out += T(x0 - 50, y0 + 22, f"No element called: {Sm['none_fraction'] * 100:.0f}% of promoters (shuffled copies {Sm['none_fraction_shuffled'] * 100:.0f}%)", 10, MUTED)
    for k, f in enumerate((0.1, 0.45, 0.8)):
        out += f'<circle cx="{lx + 110 + k * 18:.1f}" cy="{ly0 + 18:.1f}" r="6" fill="{arch_colour(f)}" fill-opacity="0.85"></circle>'
    out += T(lx + 98, ly0 + 22, "CpG share", 9.5, MUTED, "end") + T(lx + 156, ly0 + 22, "low → high", 9.5, MUTED)
    out += f'<circle cx="{lx + 110:.1f}" cy="{ly0 + 36:.1f}" r="6" fill="#ffffff" stroke="#5b574d" stroke-dasharray="3 2"></circle>' + T(lx + 122, ly0 + 40, "< 2× shuffled", 9.5, MUTED)
    out += T(lx + 98, ly0 + 40, "(×n) promoters / shuffled", 9.5, MUTED, "end")
    bottom = (ys[3] + max(n["r"] for n in by[3]) + 40) if by[3] else band2 + 10
    return out, bottom


def fig_s1():
    Sm = json.load(open("architecture_summary_C1ev_v7.json"))
    out = panel_label(20, 28, "", "promoters called per element (family-weighted): real promoters (bar) and their shuffled copies (tick)")
    bx, by0, bw = 160, 60, 500
    els = sorted([e for e in Sm["elements"] if e["weighted"] > 0], key=lambda e: -e["weighted"])
    vmax = 2400
    for i, e in enumerate(els):
        yy = by0 + i * 22
        out += T(bx - 8, yy + 11, e["element"], 10.5, INK, "end")
        out += R(bx, yy, e["weighted"] / vmax * bw, 14, BLUE, extra=' fill-opacity="0.8"')
        xs_ = bx + e["shuffled"] / vmax * bw
        out += L(xs_, yy - 2, xs_, yy + 16, ORANGE, 2.2)
        out += T(bx + e["weighted"] / vmax * bw + 6, yy + 11, f"{e['weighted']:,.0f}  ({100 * e['true_fraction']:.0f}% above chance)", 9.5, MUTED)
    for v in (0, 1000, 2000):
        xx = bx + v / vmax * bw
        out += L(xx, by0 + len(els) * 22 + 2, xx, by0 + len(els) * 22 + 6, "#8a8578") + T(xx, by0 + len(els) * 22 + 18, f"{v:,}", 9.5, MUTED, "middle")
    out += T(bx + bw / 2, by0 + len(els) * 22 + 34, "family-weighted promoters", 10, MUTED, "middle")
    body = header("S1", "S1 Fig · Element calls are mostly above chance; GC box calls are the noisiest",
                  "Each element is called only inside its positional zone; shuffled copies of the same promoters give the number of calls expected by chance.",
                  "Final weights (17 Oct) · to draw")
    h = by0 + len(els) * 22 + 50
    body += svg(h, out, "Element calls in real and shuffled promoters")
    body += footer("datasets/profiles/architecture_summary_C1ev_v7.json<br>architecture_calls_C1ev_v7.tsv<br>script: hpsad_architecture_v7.py", ["Supplementary only (decided 16 Oct)."])
    return page("S1 Fig sketch", h + 300, body)

boards = [("Main.dc.html", "Overview · story arc", main_board, 860), ("Fig1.dc.html", "Fig 1 · Anchors", fig1, 700),
          ("Fig2.dc.html", "Fig 2 · Cohorts & families", fig2, 1150), ("Fig3.dc.html", "Fig 3 · Test bench", fig3, 860),
          ("Fig4.dc.html", "Fig 4 · Detection limits", fig4, 740), ("Fig5.dc.html", "Fig 5 · C1 motifs", fig5, None),
          ("Fig6.dc.html", "Fig 6 · Start-site region", fig_tss, None),
          ("Fig7.dc.html", "Fig 7 · Across cohorts", fig6, 1030), ("Fig8.dc.html", "Fig 8 · −17 scar & TSDs", fig7, 1040),
          ("FigS1.dc.html", "S1 Fig · Element calls", fig_s1, None)]
heights = {}
for fn, title, f, h in boards:
    CUR = ("S1" if fn.startswith("FigS1") else fn[3]) if fn.startswith("Fig") else None
    html = f()
    open(f"{ROOT}/{fn}", "w").write(html)
    import re
    heights[fn] = int(re.search(r'"height":(\d+)', html).group(1))
# layout: two columns
pos, yl, yr = {}, 0, 0
for i, (fn, title, f, h) in enumerate(boards):
    hh = heights[fn]
    if i % 2 == 0:
        pos[fn] = (0, yl); yl += hh + 120
    else:
        pos[fn] = (1120, yr); yr += hh + 120
idx = {"v": 3, "createdOnFiles": {"v": 1, "at": "2026-10-12T22:45:00Z"}, "title": "HPSAD – figure storyboard",
       "launch": {"view": "canvas"}, "pages": [],
       "boards": {fn: {"x": pos[fn][0], "y": pos[fn][1], "w": W, "h": heights[fn], "title": title} for fn, title, f, h in boards},
       "order": [fn for fn, *_ in boards], "notes": {}, "designSystems": []}
json.dump(idx, open(f"{ROOT}/canvas.json", "w"), ensure_ascii=False)
print(heights)
