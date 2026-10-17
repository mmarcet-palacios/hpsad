# run_final_v7.ps1 - final reruns with corrected family weights (16 Oct 2026)
# Family size now counts distinct loci (a main and an alternative start site of one gene are the same DNA).
# Order: families -> start-site copies -> Fig 5 discovery (6 replicates) -> profiles + syntax -> Fig 6 discovery
#        -> start-site test + position-free check + architectures -> Fig 7 discovery (10 cohorts).
# Stops at the first failing step. Log: datasets\qc\run_final_v7.log
$ErrorActionPreference = 'Stop'
Set-Location C:\Promoter\HPSAD_v7
$py = 'C:\Promoter\venv_gpu\Scripts\python.exe'
$log = 'datasets\qc\run_final_v7.log'
$t0 = Get-Date

function Step($name, [string[]]$cmd) {
    $ErrorActionPreference = 'Continue'      # python warnings on stderr must not stop the run; the exit code decides
    $s = Get-Date
    "`n===== $name  ($(Get-Date -Format 'HH:mm'))" | Tee-Object -FilePath $log -Append
    & $py -u @cmd 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $log -Append
    if ($LASTEXITCODE) { throw "$name failed (exit $LASTEXITCODE) - see $log" }
    "===== $name done in $([int]((Get-Date) - $s).TotalMinutes) min" | Tee-Object -FilePath $log -Append
}

# keep the previous family files for comparison
foreach ($f in 'hpsad_families_v7.tsv', 'hpsad_families_v7_summary.json', 'hpsad_family_pairs_v7.tsv',
               'hpsad_families_v7_down100.tsv', 'hpsad_families_v7_down100_summary.json') {
    $p = "datasets\families\$f"
    if (Test-Path $p) { Copy-Item $p ($p -replace '(\.\w+)$', '_before_locusfix$1') -Force }
}

Step 'families (upstream copies)'      @('scripts\hpsad_families_v7.py')
Step 'families (start-site copies)'    @('scripts\hpsad_families_core_v7.py')
Step 'Fig 5 discovery, 6 replicates'   @('scripts\run_cohorts_parallel_v7.py', '--gpu', '--jobs', '22', '--nulls', '10', '--max-motifs', '40',
                                         '--only', 'C1_evidence', '--replicates', '_f1:61000,_f2:62000,_f3:63000,_f4:64000,_f5:65000,_f6:66000')
Step 'Fig 5 profiles'                  @('scripts\hpsad_profile_v7.py', '--tags', '_f1,_f2,_f3')
Step 'Fig 5 syntax (TATA groups)'      @('scripts\hpsad_syntax_v7.py', '--tags', '_f1,_f2,_f3', '--write-ids')
Step 'Fig 5 syntax'                    @('scripts\hpsad_syntax_v7.py', '--tags', '_f1,_f2,_f3')
Step 'Fig 6 discovery (-300..+100)'    @('scripts\run_cohorts_parallel_v7.py', '--gpu', '--jobs', '22', '--nulls', '10', '--max-motifs', '60',
                                         '--only', 'C1ev_down,C1ev_down_tata_pos,C1ev_down_tata_neg_cpg,C1ev_down_tata_neg_noncpg',
                                         '--tag', '_f', '--seed-base', '67000')
Step 'Fig 6 start-site test'           @('scripts\hpsad_downstream_v7.py', '--tag', '_f')
Step 'position-free check'             @('scripts\hpsad_kmer_enrichment_v7.py', '--tag', '_f')
Step 'Fig 7B architectures'            @('scripts\hpsad_architecture_v7.py', '--tags', '_f1,_f2,_f3', '--down-tag', '_f')
Step 'Fig 7 discovery (10 cohorts)'    @('scripts\run_cohorts_parallel_v7.py', '--gpu', '--jobs', '22', '--nulls', '5', '--max-motifs', '60',
                                         '--only', 'C1_pol2_d,C1_protein_coding_d,C1_non_coding_d,C1_small_rna_d,C1_pseudogene_d,C2_d,C2b_d,C3_d,C4_pol2_d,C4_small_rna_d',
                                         '--tag', '_f', '--seed-base', '68000')
"`nALL DONE in $([int]((Get-Date) - $t0).TotalMinutes) min" | Tee-Object -FilePath $log -Append
