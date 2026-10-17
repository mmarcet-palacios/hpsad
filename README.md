# HPSAD — Hierarchical Promoter Sequence Architecture Discovery

Code, result tables and figures for:

> Marcet-Palacios M. *De novo positional discovery of core promoter elements across 69,705 human promoters.* (2026, submitted)

HPSAD asks one question of every 6–8-letter DNA word: **does it pile up at one fixed distance from the transcription start site, more than in shuffled copies of the same promoters?** No motif database is used for discovery; known factors are named only afterwards.

Think of it as lining up thousands of promoters by their start site, like runners on a starting line, and looking for words that always stand at the same spot behind the line.

## Repository layout

```
scripts/            all analysis code (Python 3.11+, PowerShell run files for Windows)
data/               result tables used in the paper (mirrors datasets/ of a full run)
  cohorts/          anchors (start sites) and cohort summaries
  families/         pseudogene/paralog family groups and weights (1 / family size)
  motifs/           discovered elements, per cohort and per replicate
  profiles/         positional profiles, TATA syntax, start-site tests, architectures
  discovery/        word-level discovery and test-bench results (Fig 3)
  qc/               anchor checks, target-site duplications, L1 footprint
figures/
  png/              Figures 1–8 and S1 as submitted
  svg/              the same figures as editable vector files
  source/           draw_figures.py and the JSON data each panel is drawn from
```

## Inputs (not included; downloaded by the scripts or by hand)

| Input | Source | Used by |
|---|---|---|
| GRCh38.p14 genome FASTA `GCF_000001405.40_GRCh38.p14_genomic.fna.gz` | NCBI RefSeq | extract_cohort_sequences_v7.py, hpsad_tsd_v7.py |
| GRCh38.p14 annotation `GCF_000001405.40_GRCh38.p14_genomic.gff.gz` | NCBI RefSeq | build_cohorts_v7.py, hpsad_tsd_v7.py |
| EPDnew human promoters v006 | https://epd.expasy.org/ (download_epdnew.py) | compare_tss_anchors.py, build_cohorts_v7.py |
| FANTOM5 CAGE peaks, hg38 | https://fantom.gsc.riken.jp/5/ (download_fantom5_cage.py) | match_cage_anchors_v7.py |

Place the two NCBI files in `datasets/raw/`. The large intermediate sequence tables (`hpsad_segments_v7.tsv.gz`, `hpsad_segments_v7_down100.tsv.gz`) are archived in the Zenodo record linked below.

## Run order

All scripts take `--root` (default: current folder) and read/write under `datasets/`.

1. **Anchors and cohorts** (Fig 1, Fig 2A)
   `download_epdnew.py` → `download_fantom5_cage.py` → `compare_tss_anchors.py` → `match_cage_anchors_v7.py` → `build_cohorts_v7.py --cage` → `extract_cohort_sequences_v7.py` → `extract_downstream_v7.py --down 100`
2. **Families and weights** (Fig 2B–D) — `hpsad_families_v7.py`, `hpsad_families_core_v7.py`; LINE-1 checks (Fig 8): `hpsad_tsd_v7.py`, `hpsad_l1_footprint_v7.py`
3. **Discovery, profiles, start-site test, architectures** (Figs 5–7) — `run_final_v7.ps1`, then `run_followup_f.ps1`
4. **Test benches** — Fig 3: `hpsad_discovery_v7.py --bench` (and the hard bench in `run_followup_f.ps1`); Fig 4: `run_bench_f.ps1`
5. **Figures** — `cd figures/source; python draw_figures.py` writes one HTML page per figure to `out/project/`; `node render_png.js <page.html> <out.png> …` (Playwright) renders each figure to PNG at 3×

`hpsad_coverage_v7.py` is an optional check of how many promoters the elements explain. `hpsad_motifs_v7.py` is the shared engine (word scan, shuffled nulls, family weighting, motif building). `hpsad_gpu_check.py` confirms the optional GPU path gives the same numbers as the CPU path.

## Figure → script → data

| Figure | Scripts | Data |
|---|---|---|
| 1 Annotated vs experimental start sites | compare_tss_anchors.py | data/qc/tss_anchor_comparison* |
| 2 Cohorts and families | build_cohorts_v7.py, hpsad_families_v7.py, hpsad_families_core_v7.py | data/cohorts/*summary.json, data/families |
| 3 Word-level discovery and bench | hpsad_discovery_v7.py | data/discovery |
| 4 Motif test bench | hpsad_motifs_v7.py --bench (run_bench_f.ps1) | data/motifs/hpsad_motifs_bench_zones_f_summary.json |
| 5 Elements upstream of the start site | run_cohorts_parallel_v7.py, hpsad_stability_v7.py, hpsad_profile_v7.py, hpsad_syntax_v7.py | data/motifs/*C1_evidence_f1–f6*, data/profiles/profile_*_6rep*, syntax_*_6rep* |
| 6 Elements at the start site | hpsad_downstream_v7.py, hpsad_consensus_profile_v7.py, hpsad_kmer_enrichment_v7.py | data/motifs/*C1ev_down*, data/profiles/downstream_*, kmer_enrichment_* |
| 7 Architectures and cohorts | hpsad_architecture_v7.py, hpsad_motif_annotation_v7.py, figures/source/assign_fig7_elements.py | data/profiles/architecture_*, data/motifs/*_d_f*, figures/source/fig7_f.json |
| 8 The motif at −17 and target-site duplications | hpsad_l1_footprint_v7.py, hpsad_tsd_v7.py | data/qc/l1_footprint_v7.tsv, tsd_v7_* |
| S1 Element calls, real vs shuffled | hpsad_architecture_v7.py | data/profiles/architecture_calls_C1ev_v7.tsv, architecture_summary_C1ev_v7.json |

## Hardware

The paper's runs used Windows 11, 24 CPU threads and an NVIDIA GPU (CuPy). The CPU path gives identical results; the full discovery run takes roughly a day on CPU.

## Licence and citation

Code: MIT (see `LICENSE`). Data tables: CC BY 4.0. Please cite the article and the archived release (`CITATION.cff`).

Archive DOI: [10.5281/zenodo.23441932](https://doi.org/10.5281/zenodo.23441932) (includes the promoter sequence tables)
