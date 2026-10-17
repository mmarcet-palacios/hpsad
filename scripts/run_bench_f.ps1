# run_bench_f.ps1 - Fig 4 detection-limit bench on the final (locus-based) family weights (17 Oct 2026)
# Old null-motif caches for C1 were made with the old weights; they are moved aside so the nulls are rerun.
# Output: datasets\motifs\hpsad_motifs_bench_zones_f_summary.json, hpsad_motifs_C1_planted_zones_f.tsv
# Log: datasets\qc\run_bench_f.log
$ErrorActionPreference = 'Stop'
Set-Location C:\Promoter\HPSAD_v7
$py = 'C:\Promoter\venv_gpu\Scripts\python.exe'
$log = 'datasets\qc\run_bench_f.log'
$t0 = Get-Date
$old = 'datasets\motifs\nulls\old_weights'
New-Item -ItemType Directory -Force $old | Out-Null
Get-ChildItem 'datasets\motifs\nulls' -Filter 'motifs_null_C1_*.json' | Move-Item -Destination $old -Force
"moved old C1 null-motif caches to $old" | Tee-Object -FilePath $log
$env:HPSAD_GPU = '1'
$ErrorActionPreference = 'Continue'
& $py -u scripts\hpsad_motifs_v7.py --bench --cohort C1 --nulls 5 --jobs 18 --tag _f 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $log -Append
if ($LASTEXITCODE) { throw "bench failed (exit $LASTEXITCODE) - see $log" }
"`nALL DONE in $([int]((Get-Date) - $t0).TotalMinutes) min" | Tee-Object -FilePath $log -Append
