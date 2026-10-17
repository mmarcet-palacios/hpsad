# run_followup_f.ps1 - steps after run_final_v7.ps1 that produce the Fig 5-7 and S1 tables used in the paper
# (six-replicate elements, start-site profiles with the YY1 matrix of the main start-site run, consensus profiles,
#  and the word-level test bench of Fig 3). Stops at the first failing step.
$ErrorActionPreference = 'Stop'
Set-Location C:\Promoter\HPSAD_v7
$py = 'C:\Promoter\venv_gpu\Scripts\python.exe'
function Step($cmd) { $ErrorActionPreference = 'Continue'; & $py -u @cmd; if ($LASTEXITCODE) { throw "failed: $cmd" } }
Step @('scripts\hpsad_profile_v7.py', '--tags', '_f1,_f2,_f3,_f4,_f5,_f6', '--min-reps', '3', '--suffix', '_6rep')
Step @('scripts\hpsad_syntax_v7.py', '--tags', '_f1,_f2,_f3,_f4,_f5,_f6', '--min-reps', '3', '--suffix', '_6rep')
Step @('scripts\hpsad_downstream_v7.py', '--tag', '_f')
Step @('scripts\hpsad_consensus_profile_v7.py', '--tag', '_f')
Step @('scripts\hpsad_architecture_v7.py', '--tags', '_f1,_f2,_f3,_f4,_f5,_f6', '--min-reps', '3', '--down-tag', '_f')
Step @('scripts\hpsad_discovery_v7.py', '--bench')
Step @('scripts\hpsad_discovery_v7.py', '--bench', '--seed', '23', '--fidelity', '0.75', '--fractions', '0.0005,0.001,0.002,0.005', '--jitters', '0,15,30', '--tag', '_hard')
"done"
