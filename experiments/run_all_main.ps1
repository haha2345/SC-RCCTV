# Full re-run of every experiment of the paper on ONE machine (the main machine), so that all results share the same
# hardware, MATLAB release, BLAS and inputs. Outputs go to a fresh folder (default results_main\); results\ is only
# read (noisy inputs) and never modified, so nothing of the earlier runs is lost.
#
#   powershell -ExecutionPolicy Bypass -File run_all_main.ps1 -Matlab "C:\Program Files\MATLAB\R2022b\bin\matlab.exe"
#
# Phases (default: all, in this order; every phase is resumable, so the script can be stopped and restarted):
#   check   MATLAB release, toolboxes, third-party code, Python; SHA-256 of all inputs against inputs_manifest.csv
#   s3ttv   S3TTV (seed 0, cases 1-8). Default: keep the S3TTV runs that were already made on this machine
#           (results\sim_s3*.csv are copied); -WithS3TTV re-runs them (6-10 h on CPU)
#   tune    t-CTV sparse-weight multiplier on the development cubes; the best multiplier is used for all t-CTV runs
#   sim     all methods on the 116 simulated cubes (t-CTV with the tuned multiplier on all six realizations)
#   real    Indian Pines and Urban, all methods
#   extra   rank sensitivity (Indian Pines) and the Python checks of the stripe false-alarm rate and noise-level settling
#   idle    wall-clock timing and full-size scenes; ONE MATLAB at a time, the machine must not be used meanwhile
#   analyze real-data metrics, SVM classification, all figures
#   summary completeness report (results_main\_SUMMARY.txt)
# Example: -Phase sim runs one phase only; -Par 6 sets the number of parallel MATLAB processes outside the idle phase.
param(
  [string]$Matlab = 'matlab',
  [string]$Python = 'python',
  [ValidateSet('all', 'check', 's3ttv', 'tune', 'sim', 'real', 'extra', 'idle', 'analyze', 'summary')][string]$Phase = 'all',
  [int]$Par = 6,
  [string]$ResultsDir = 'results_main',
  [switch]$WithS3TTV,
  [ValidateSet('auto', '0', '1')][string]$Gpu = 'auto',
  [switch]$Smoke                      # tiny subset, only to test the plumbing (use with a scratch -ResultsDir)
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path        # ...\new_method_2026-10-12
$mdir = Join-Path $root 'matlab'
$out = if ([IO.Path]::IsPathRooted($ResultsDir)) { $ResultsDir } else { Join-Path $root $ResultsDir }
$logd = Join-Path $out 'logs'
foreach ($d in $out, $logd, (Join-Path $out 'real'), (Join-Path $out 'figs')) { New-Item -ItemType Directory -Force $d | Out-Null }
$env:RESULTS_DIR = $out                                         # read by matlab\res_dirs.m: every output goes here
$cores = [Environment]::ProcessorCount
$thr = [math]::Max(2, [int][math]::Floor($cores / [math]::Max(1, $Par)))
$runlog = Join-Path $out '_run_log.txt'

function Log([string]$msg) { $l = "[$(Get-Date -Format 'MM-dd HH:mm:ss')] $msg"; Write-Host $l; Add-Content -Path $runlog -Value $l }

$JOBKEYS = 'SIM_CUBE', 'SIM_METHS', 'SIM_SEEDS', 'SIM_CASES', 'SIM_TAG', 'TCTV_MULT', 'REAL_DS', 'REAL_METH', 'TUNE_CUBES',
           'TUNE_MULTS', 'TUNE_TAG', 'LARGE_JOBS', 'LARGE_METHS', 'S3TTV_MAXIT', 'S3TTV_GPU', 'SKIP_CLS'
function Start-MlJob($job) {
  foreach ($k in $JOBKEYS) { Remove-Item "env:$k" -ErrorAction SilentlyContinue }
  foreach ($kv in $job.Env.GetEnumerator()) { Set-Item -Path "env:$($kv.Key)" -Value $kv.Value }
  $pre = if ($job.Threads -gt 0) { "maxNumCompThreads($($job.Threads)); " } else { '' }
  $stamp = Get-Date -Format 'MMdd_HHmmss'
  $p = Start-Process -FilePath $Matlab -ArgumentList "-batch `"$pre$($job.Cmd)`"" -WorkingDirectory $mdir `
        -RedirectStandardOutput (Join-Path $logd "$($job.Name)_$stamp.log") -RedirectStandardError (Join-Path $logd "$($job.Name)_$stamp.err") `
        -PassThru -WindowStyle Hidden
  $null = $p.Handle                       # keep the handle, otherwise ExitCode is empty after the exit
  Log "start $($job.Name) (pid $($p.Id))"
  return $p
}
function Invoke-Pool($jobs, [int]$par) {
  $run = @()
  foreach ($j in @($jobs)) {
    while (@($run | Where-Object { -not $_.P.HasExited }).Count -ge $par) { Start-Sleep -Seconds 10 }
    $run += [pscustomobject]@{ P = (Start-MlJob $j); Name = $j.Name }
  }
  foreach ($r in $run) { $r.P.WaitForExit(); Log "end   $($r.Name) (exit code $($r.P.ExitCode))" }
  foreach ($k in $JOBKEYS) { Remove-Item "env:$k" -ErrorAction SilentlyContinue }
}
function Invoke-One($name, $cmd, $envs = @{}, [int]$threads = 0) {
  Invoke-Pool @(@{ Name = $name; Cmd = $cmd; Env = $envs; Threads = $threads }) 1
}
function Get-PythonOk {
  if (-not (Get-Command $Python -ErrorAction SilentlyContinue)) { return $false }
  $o = Join-Path $logd 'python_check.txt'
  $p = Start-Process -FilePath $Python -ArgumentList '-c "import numpy, scipy, sys; print(sys.version); print(numpy.__version__, scipy.__version__)"' `
        -RedirectStandardOutput $o -RedirectStandardError "$o.err" -PassThru -WindowStyle Hidden -WorkingDirectory $root
  $null = $p.Handle
  $p.WaitForExit()
  return ($p.ExitCode -eq 0)
}

$SEEDS = '0,1,2,3,4,5'; $CASES = '1,2,3,4,5,6,7,8,9,10'; $FAST = 'RCCTV,OLD,OLDband,SC,SCv,SC_scalar,SC_noT,SC_noG,SC_fix'
$CUBES = @(@('pavia120', 'P'), @('belt100', 'B'))
if ($Smoke) { $SEEDS = '0'; $CASES = '8'; $FAST = 'RCCTV,SC,SC_fix'; $CUBES = @(, @('pavia120', 'P')) }

# ---------------------------------------------------------------------------------------------------- check
function Phase-Check {
  Log "check: results go to $out; $cores logical processors, $Par parallel jobs x $thr threads"
  Invoke-One 'env_check' 'env_check'
  $envf = Join-Path $out 'env_main.txt'
  if (-not (Test-Path $envf) -or -not (Select-String -Path $envf -Pattern 'STATUS OK' -Quiet)) { throw "MATLAB environment check failed, see $envf and $logd" }
  $cpu = Get-CimInstance Win32_Processor | Select-Object -First 1
  $os = Get-CimInstance Win32_OperatingSystem
  Add-Content $envf ("CPU {0}; {1} cores / {2} logical; RAM {3:N1} GB; {4} {5}" -f $cpu.Name, $cpu.NumberOfCores, $cpu.NumberOfLogicalProcessors, ($os.TotalVisibleMemorySize / 1MB), $os.Caption, $os.Version)
  if (Get-PythonOk) { Add-Content $envf ("Python " + ((Get-Content (Join-Path $logd 'python_check.txt')) -join ' | ')) }
  else { Add-Content $envf 'Python NOT AVAILABLE (numpy, scipy): the false-alarm and settling checks will be skipped'; Log 'WARNING: Python with numpy and scipy not found' }
  $man = Import-Csv (Join-Path $root 'inputs_manifest.csv')
  $bad = @(); $warn = @()
  foreach ($m in $man) {
    $p = Join-Path $root $m.path
    $msg = if (-not (Test-Path $p)) { "missing $($m.path)" } elseif ((Get-FileHash $p -Algorithm SHA256).Hash -ne $m.sha256) { "differs $($m.path)" } else { $null }
    if ($msg) { if ($m.kind -eq 'abort') { $bad += $msg } else { $warn += $msg } }
  }
  $warn | ForEach-Object { Log "WARNING input $_ (large-scale input; the run continues with this file)" }
  if ($bad) { $bad | Select-Object -First 20 | ForEach-Object { Log "ERROR input $_" }; throw "$($bad.Count) input files are missing or differ from inputs_manifest.csv" }
  Log "check: $($man.Count) input files verified (SHA-256)"
}

# ---------------------------------------------------------------------------------------------------- S3TTV
function Phase-S3 {
  if (-not $WithS3TTV) {
    foreach ($t in 'sim_s3P.csv', 'sim_s3B.csv') {
      $dst = Join-Path $out $t; $src = Join-Path $root "results\$t"
      if (Test-Path $dst) { continue }
      if (Test-Path $src) { Copy-Item $src $dst; Log "s3ttv: kept the S3TTV runs already made on this machine ($t copied from results\)" }
      else { Log "WARNING: results\$t not found; run again with -WithS3TTV" }
    }
    return
  }
  $useGpu = $Gpu
  if ($Gpu -eq 'auto') {
    $o = Join-Path $logd 'gpu_probe.txt'
    $p = Start-Process -FilePath $Matlab -ArgumentList "-batch `"try, disp(gpuDeviceCount('available')>0), catch, disp(0), end`"" -RedirectStandardOutput $o -PassThru -WindowStyle Hidden -WorkingDirectory $mdir
    $p.WaitForExit(); $useGpu = ((Get-Content $o | Where-Object { $_.Trim() -in '0', '1' } | Select-Object -Last 1) -as [string]).Trim()
    if ($useGpu -notin '0', '1') { $useGpu = '0' }
  }
  $maxit = if ($useGpu -eq '1') { '20000' } else { '1000' }
  Log "s3ttv: re-run, gpu=$useGpu, maxiter=$maxit"
  $jobs = foreach ($cb in $CUBES) { @{ Name = "s3$($cb[1])"; Cmd = 'sim_run_new'; Threads = $thr
      Env = @{ SIM_CUBE = $cb[0]; SIM_METHS = 'S3TTV'; SIM_SEEDS = '0'; SIM_CASES = '8,5,6,7,4,1,2,3'; SIM_TAG = "s3$($cb[1])"; S3TTV_GPU = $useGpu; S3TTV_MAXIT = $maxit } } }
  Invoke-Pool $jobs ([math]::Min($Par, 2))
}

# ---------------------------------------------------------------------------------------------------- t-CTV tuning
function Get-BestMult {
  if ($Smoke) { return '0.5' }
  $fs = @(Get-ChildItem $out -Filter 'tune_tctv*.csv')
  if (-not $fs) { throw 'tune_tctv*.csv missing: run -Phase tune first' }
  $rows = @($fs | ForEach-Object { Import-Csv $_.FullName })
  $nm = @($rows | Select-Object -ExpandProperty mult -Unique).Count
  # paired comparison on the (cube, case) cells that were run with every multiplier (same rule as make_rsl.py)
  $cells = @($rows | Group-Object cube, case | Where-Object { @($_.Group | Select-Object -ExpandProperty mult -Unique).Count -eq $nm } | ForEach-Object { $_.Name })
  $sel = $rows | Where-Object { $cells -contains "$($_.cube), $($_.case)" }
  $means = $sel | Group-Object mult | ForEach-Object { [pscustomobject]@{ mult = $_.Name; mean = ($_.Group | ForEach-Object { [double]$_.mpsnr } | Measure-Object -Average).Average; n = $_.Count } }
  $best = $means | Sort-Object mean -Descending | Select-Object -First 1
  $means | Sort-Object { [double]$_.mult } | ForEach-Object { Log ("tune: multiplier {0,-6} mean MPSNR {1:N3} dB over {2} cells" -f $_.mult, $_.mean, $_.n) }
  if (@($cells).Count -lt 30) { Log "WARNING: only $(@($cells).Count) of 30 development cells are complete" }
  Set-Content (Join-Path $out 'tctv_mult.txt') $best.mult
  return $best.mult
}
function Phase-Tune {
  if ($Smoke) { Log 'tune: skipped in smoke mode (multiplier 0.5)'; return }
  $mults = '0.125,0.25,0.5,1,2,4'
  $jobs = foreach ($c in 'wdc1', 'wdc2', 'indian') { @{ Name = "tune_$c"; Cmd = 'tune_tctv'; Env = @{ TUNE_CUBES = $c; TUNE_MULTS = $mults; TUNE_TAG = $c }; Threads = $thr } }
  Invoke-Pool $jobs $Par
  Log "tune: best multiplier $(Get-BestMult)"
}

# ---------------------------------------------------------------------------------------------------- simulated data
function Phase-Sim([string]$mult) {
  $jobs = @()
  foreach ($cb in $CUBES) {
    $c = $cb[0]; $s = $cb[1]
    if (-not $Smoke) {
      $jobs += @{ Name = "sim_tcd${s}_a"; Cmd = 'sim_run_new'; Threads = $thr; Env = @{ SIM_CUBE = $c; SIM_METHS = 'tCTVd'; SIM_SEEDS = $SEEDS; SIM_CASES = '1,2,3,4,5'; SIM_TAG = "tcd${s}_a"; TCTV_MULT = $mult } }
      $jobs += @{ Name = "sim_tcd${s}_b"; Cmd = 'sim_run_new'; Threads = $thr; Env = @{ SIM_CUBE = $c; SIM_METHS = 'tCTVd'; SIM_SEEDS = $SEEDS; SIM_CASES = '6,7,8,9,10'; SIM_TAG = "tcd${s}_b"; TCTV_MULT = $mult } }
      $jobs += @{ Name = "sim_lr$s"; Cmd = 'sim_run'; Threads = $thr; Env = @{ SIM_CUBE = $c; SIM_METHS = 'LRTDTV'; SIM_SEEDS = $SEEDS; SIM_CASES = $CASES; SIM_TAG = "lr$s" } }
      $jobs += @{ Name = "sim_tc$s"; Cmd = 'sim_run'; Threads = $thr; Env = @{ SIM_CUBE = $c; SIM_METHS = 'tCTV'; SIM_SEEDS = '0'; SIM_CASES = '1,2,3,4,5,6,7,8'; SIM_TAG = "tc$s" } }
      $jobs += @{ Name = "sim_e3$s"; Cmd = 'sim_run'; Threads = $thr; Env = @{ SIM_CUBE = $c; SIM_METHS = 'E3DTV'; SIM_SEEDS = $SEEDS; SIM_CASES = $CASES; SIM_TAG = "e3$s" } }
    }
    $jobs += @{ Name = "sim_fh$s"; Cmd = 'sim_run_new'; Threads = $thr; Env = @{ SIM_CUBE = $c; SIM_METHS = 'FastHyMix,L1HyMixDe'; SIM_SEEDS = $SEEDS; SIM_CASES = $CASES; SIM_TAG = "fh$s" } }
    $jobs += @{ Name = "sim_fast$s"; Cmd = 'sim_run'; Threads = $thr; Env = @{ SIM_CUBE = $c; SIM_METHS = $FAST; SIM_SEEDS = $SEEDS; SIM_CASES = $CASES; SIM_TAG = "fast$s" } }
  }
  Log "sim: $($jobs.Count) jobs (t-CTV multiplier $mult)"
  Invoke-Pool $jobs $Par
}

# ---------------------------------------------------------------------------------------------------- real data
function Phase-Real([string]$mult) {
  $jobs = @()
  $dss = if ($Smoke) { @('indian220') } else { @('indian220', 'urban200') }
  $m1 = if ($Smoke) { @('SC') } else { @('RCCTV', 'SC', 'SC_noT', 'SC_noG', 'SC_scalar', 'SC_fix', 'E3DTV', 'LRTDTV') }
  $m2 = if ($Smoke) { @('FastHyMix') } else { @('FastHyMix', 'L1HyMixDe') }
  foreach ($ds in $dss) {
    foreach ($m in $m1) { $jobs += @{ Name = "real_${ds}_$m"; Cmd = 'real_run'; Threads = $thr; Env = @{ REAL_DS = $ds; REAL_METH = $m } } }
    foreach ($m in $m2) { $jobs += @{ Name = "real_${ds}_$m"; Cmd = 'real_run_new'; Threads = $thr; Env = @{ REAL_DS = $ds; REAL_METH = $m; TCTV_MULT = $mult } } }
  }
  if (-not $Smoke) {
    $jobs += @{ Name = 'real_indian220_tCTV'; Cmd = 'real_run'; Threads = $thr; Env = @{ REAL_DS = 'indian220'; REAL_METH = 'tCTV' } }
    $jobs += @{ Name = 'real_indian220_tCTVd'; Cmd = 'real_run_new'; Threads = $thr; Env = @{ REAL_DS = 'indian220'; REAL_METH = 'tCTVd'; TCTV_MULT = $mult } }
  }
  $jobs = @($jobs | Where-Object { -not (Test-Path (Join-Path $out ("real\" + $_.Env.REAL_DS + '_' + $_.Env.REAL_METH + '.mat'))) })
  Log "real: $($jobs.Count) jobs to run"
  if ($jobs.Count) { Invoke-Pool $jobs $Par }
}

# ---------------------------------------------------------------------------------------------------- extra checks
function Phase-Extra {
  $rs = Join-Path $out 'rank_sens.csv'
  if ($Smoke) { Log 'extra: rank sensitivity skipped in smoke mode' }
  elseif (-not ((Test-Path $rs) -and (@(Get-Content $rs).Count -ge 7))) { Invoke-One 'rank_sens' 'rank_sens' @{} $thr }
  if (-not (Get-PythonOk)) { Log 'WARNING: Python not available, check_fa / check_fp / theory_checks skipped'; return }
  $scripts = if ($Smoke) { @('check_fa') } else { @('check_fa', 'check_fp', 'theory_checks') }
  foreach ($py in $scripts) {
    $o = Join-Path $out "$py.txt"
    if ((Test-Path $o) -and (Get-Item $o).Length -gt 0) { continue }
    Log "extra: python $py.py"
    $p = Start-Process -FilePath $Python -ArgumentList "`"$(Join-Path $root "$py.py")`"" -WorkingDirectory $root `
          -RedirectStandardOutput $o -RedirectStandardError (Join-Path $logd "$py.err") -PassThru -WindowStyle Hidden
    $null = $p.Handle
    $p.WaitForExit()
    if ($p.ExitCode -ne 0) { Log "WARNING: $py.py failed (exit $($p.ExitCode)), see $logd\$py.err"; Remove-Item $o -ErrorAction SilentlyContinue }
  }
}

# ---------------------------------------------------------------------------------------------------- idle-machine timing
function Phase-Idle([string]$mult) {
  if ($Smoke) { Log 'idle: skipped in smoke mode'; return }
  $other = @(Get-Process MATLAB -ErrorAction SilentlyContinue)
  if ($other.Count) { Log "WARNING: $($other.Count) MATLAB process(es) already running; the timing is only valid on an otherwise idle machine" }
  Log 'idle: timing and full-size scenes, one MATLAB at a time with default threads (do not use the machine)'
  Invoke-One 'timing' 'timing'
  Invoke-One 'timing_new' 'timing_new' @{ TCTV_MULT = $mult }
  foreach ($m in 'RCCTV', 'SC', 'FastHyMix', 'L1HyMixDe', 'E3DTV', 'LRTDTV', 'tCTVd') {     # one process per method
    Invoke-One "large_$m" 'large_run' @{ LARGE_JOBS = 'paviaU_c8,paviaU_c10,urban307'; LARGE_METHS = $m; TCTV_MULT = $mult }
  }
}

# ---------------------------------------------------------------------------------------------------- analysis + figures
function Phase-Analyze {
  if ($Smoke) { Log 'analyze: skipped in smoke mode'; return }
  Invoke-One 'real_analyze' 'real_analyze'          # metrics, SVM on Indian Pines, real-data figures, sigma curve
  Invoke-One 'sim_fig' 'sim_fig'
  Invoke-One 'pareto_fig' 'pareto_fig'
}

# ---------------------------------------------------------------------------------------------------- summary
function Phase-Summary {
  $lines = @("Summary of $out  ($(Get-Date -Format 'yyyy-MM-dd HH:mm'))", '')
  $sim = @(Get-ChildItem $out -Filter 'sim_*.csv' | ForEach-Object { Import-Csv $_.FullName })
  $expect = @{ tCTV = 8; S3TTV = 8 }
  $ok = $true
  foreach ($g in ($sim | Group-Object cube, meth | Sort-Object Name)) {
    $cube, $meth = $g.Name -split ', '
    $exp = if ($expect.ContainsKey($meth)) { $expect[$meth] } elseif ($Smoke) { 1 } else { 58 }
    $flag = if ($g.Count -ge $exp) { 'ok' } elseif ($meth -eq 'FastHyMix') { 'released code crashed in some runs (expected)' } else { $ok = $false; 'INCOMPLETE' }
    $lines += ('{0,-9} {1,-10} {2,3}/{3,-3} {4}' -f $cube, $meth, $g.Count, $exp, $flag)
  }
  $need = 'tune_tctv_wdc1.csv', 'tune_tctv_wdc2.csv', 'tune_tctv_indian.csv', 'tctv_mult.txt', 'real_metrics.csv', 'real_indian_cls.csv', 'sigma_stats.csv', 'real_picks.csv', 'rank_sens.csv',
          'timing.csv', 'timing_new.csv', 'large.csv', 'check_fa.txt', 'check_fp.txt', 'env_main.txt',
          'figs\fig_sim_visual.pdf', 'figs\fig_real_visual.pdf', 'figs\fig_sigma_curve.pdf', 'figs\fig_pareto.pdf'
  $lines += ''
  foreach ($f in $need) { $e = Test-Path (Join-Path $out $f); if (-not $e) { $ok = $false }; $lines += ('{0,-26} {1}' -f $f, $(if ($e) { 'ok' } else { 'MISSING' })) }
  $errs = @(Get-ChildItem $logd -Filter '*.log' | Select-String -Pattern '^ERROR' | ForEach-Object { "$($_.Filename): $($_.Line)" })
  $lines += '', "ERROR lines in the MATLAB logs: $($errs.Count)"
  $lines += $errs | Select-Object -First 50
  $lines += '', $(if ($ok) { 'STATUS COMPLETE' } else { 'STATUS INCOMPLETE (re-run the script; finished cells are skipped)' })
  Set-Content (Join-Path $out '_SUMMARY.txt') $lines
  $lines | ForEach-Object { Write-Host $_ }
  if ($ok) { Set-Content (Join-Path $out '_COMPLETE.txt') (Get-Date -Format 'yyyy-MM-dd HH:mm') }
}

# ---------------------------------------------------------------------------------------------------- main
Log "run_all_main.ps1 -Phase $Phase -Par $Par -ResultsDir $out$(if ($WithS3TTV) { ' -WithS3TTV' })$(if ($Smoke) { ' -Smoke' })"
if ($Phase -in 'all', 'check') { Phase-Check }
if ($Phase -in 'all', 's3ttv') { Phase-S3 }
if ($Phase -in 'all', 'tune') { Phase-Tune }
$mult = if ($Phase -in 'all', 'sim', 'real', 'idle') { Get-BestMult } else { $null }
if ($Phase -in 'all', 'sim') { Phase-Sim $mult }
if ($Phase -in 'all', 'real') { Phase-Real $mult }
if ($Phase -in 'all', 'extra') { Phase-Extra }
if ($Phase -in 'all', 'idle') { Phase-Idle $mult }
if ($Phase -in 'all', 'analyze') { Phase-Analyze }
if ($Phase -in 'all', 'summary') { Phase-Summary }
Log "finished -Phase $Phase. Copy back the whole folder $out (results_main\real\*.mat is optional) to the laptop."
