# Re-measure the S3TTV run time on the CPU (as in its evaluation runs) after the full re-run of 2026-10-14, where
# timing_new.m let S3TTV pick the GPU and fell back to the slow per-block prox (75 s/iter instead of ~3.6 s/iter).
# Removes the S3TTV rows from <ResultsDir>\timing_new.csv, runs timing_new.m (only the missing rows are measured) and
# redraws the Pareto figure. The machine must be otherwise idle (a few minutes).
#   powershell -ExecutionPolicy Bypass -File fix_s3ttv_timing.ps1 -Matlab "C:\Program Files\MATLAB\R2022a\bin\matlab.exe"
param([string]$Matlab = 'matlab', [string]$ResultsDir = 'results_main')
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path        # ...\new_method_2026-10-12
$out = if ([IO.Path]::IsPathRooted($ResultsDir)) { $ResultsDir } else { Join-Path $root $ResultsDir }
$f = Join-Path $out 'timing_new.csv'
if (-not (Test-Path $f)) { throw "$f not found" }
Copy-Item $f "$f.bak_gpu" -Force                                 # keep the invalid measurement for the record
$keep = Get-Content $f | Where-Object { $_ -notmatch ',S3TTV,' }
Set-Content $f $keep
$env:RESULTS_DIR = $out
$mf = Join-Path $out 'tctv_mult.txt'
$env:TCTV_MULT = if (Test-Path $mf) { (Get-Content $mf | Select-Object -First 1).Trim() } else { '0.5' }
$logd = Join-Path $out 'logs'
foreach ($cmd in 'timing_new', 'pareto_fig') {
  $stamp = Get-Date -Format 'MMdd_HHmmss'
  $p = Start-Process -FilePath $Matlab -ArgumentList "-batch `"$cmd`"" -WorkingDirectory (Join-Path $root 'matlab') `
        -RedirectStandardOutput (Join-Path $logd "${cmd}_fix_$stamp.log") -RedirectStandardError (Join-Path $logd "${cmd}_fix_$stamp.err") -PassThru -WindowStyle Hidden
  $null = $p.Handle; $p.WaitForExit()
  Write-Host "$cmd finished (exit code $($p.ExitCode))"
}
Get-Content $f | Where-Object { $_ -match ',S3TTV,' }
'Copy back results_main\timing_new.csv and results_main\figs\fig_pareto.* (and pareto_points.mat) to the laptop.'
