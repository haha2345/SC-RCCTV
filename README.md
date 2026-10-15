# SC-RCCTV: self-calibrated subspace robust PCA for hyperspectral mixed noise and line stripes

MATLAB implementation of **SC-RCCTV**, a self-calibrated representative-coefficient correlated total variation
(RCCTV) model for removing Gaussian noise, impulse noise and column/row stripes from hyperspectral images (HSIs)
without any image-specific parameter, together with every script used to produce the experimental results of the
accompanying manuscript.

- **Model**: RCCTV (low-rank coefficient images whose spatial gradients are penalized by nuclear norms) plus a
  band-wise dense term, a sparse term and line-constant column/row stripe terms.
- **Self-calibration**: the noise level of every band is initialized by a stripe-invariant median absolute deviation
  (MAD) estimate and refined online from the residual, which turns the sparse/dense separation into a band-adaptive
  Huber loss. The stripe thresholds follow from a significance test with false-alarm probability 2Φ(−c), 0.27% for
  c = 3. No noise level or stripe radius has to be specified.
- **Solver**: ADMM with closed-form subproblems; the per-iteration cost is that of RCCTV (about 1 s for a
  120×120×99 cube).

## Quick start

```matlab
demo_sc_rcctv          % synthetic cube with band-dependent noise, impulses and stripes; no data download needed
```

Expected output (MATLAB R2022a): `MPSNR: noisy 16.20 dB -> SC-RCCTV 38.54 dB`.

## Usage

```matlab
addpath experiments/matlab
[X, info] = sc_rcctv(Y);            % Y: M x N x p noisy HSI, each band scaled to [0, 1]
[X, info] = sc_rcctv(Y, opts);      % optional settings, see below
```

`sc_rcctv.m` is self-contained and only needs `psf2otf` from the Image Processing Toolbox.

| option | default | meaning |
|---|---|---|
| `rk` | 8 | rank of the spectral subspace |
| `lam` | 6/sqrt(MN) | weight of the sparse term (as in RCCTV) |
| `cT` | 3 | significance level of the stripe test (false-alarm probability 2Φ(−c)) |
| `stripe_dirs` | `'vh'` | stripe orientations: `'v'` column-constant, `'h'` row-constant, `'vh'` both |
| `sig_online`, `online_start` | `true`, 5 | refine the noise levels from the residual from this iteration on |
| `tol`, `maxIter` | 1e-5, 100 | stopping rule |
| `scalar_sigma`, `T_mode`, `G_mode` | | ablation switches (one global σ, no stripe term `'none'`, no dense term `'none'`) |

`info` returns the iterations (`iters`), the final noise levels (`sig_k`), and the stripe (`T`), dense (`G`) and
sparse (`E`) components.

## Reproducing the experiments

```
SC-RCCTV/
  demo_sc_rcctv.m
  experiments/          scripts (MATLAB, Python, PowerShell), input checksums, results of the full run
    matlab/             sc_rcctv.m and all experiment drivers
    results_main/       every table and figure of the paper, recomputed in one run on one machine
    external/           third-party methods (not included, see external/README.md)
  code/                 third-party code (not included, see code/README.md)
  data/                 data sets (not included, see data/README.md)
  tools/                make_tctv_trpca_lam.py
```

1. Put the third-party codes and the data sets in place (`code/README.md`, `experiments/external/README.md`,
   `data/README.md`). They are not redistributed here because their licenses do not allow it or are unspecified.
2. Create the t-CTV variant used in the comparison (the official code plus one line that scales its sparse weight):
   `python tools/make_tctv_trpca_lam.py code/official_tctv/TCTV_TRPCA.m`.
3. Generate the noisy inputs (Python with numpy and scipy), from `experiments/`:
   `python gen_inputs.py`, `python gen_dev_inputs.py`, `python gen_large_inputs.py`.
   Seed 0 of cases 1–8 are the official RCCTV test inputs; seeds 1–5 and cases 9–10 are generated.
   `inputs_manifest.csv` lists the SHA-256 of every input file used for the reported results.
4. Run everything (Windows PowerShell; about 5 h on an Intel i7-14700KF):
   ```
   powershell -ExecutionPolicy Bypass -File run_all_main.ps1 -Matlab "<path to matlab.exe>"
   ```
   The script checks the environment and the input checksums, tunes the t-CTV weight on development cubes, runs all
   methods on the simulated and real data, measures the run times one MATLAB at a time, and draws the figures.
   It is resumable; `-Phase <name>` runs a single phase and `-Par <n>` sets the number of parallel MATLAB processes.
   On a fresh setup add `-WithS3TTV` (S3TTV takes about 40 min per cube on the CPU); without it, earlier S3TTV
   results in `experiments/results/` are reused.
   Two control runs added in the revision are not part of the script; they are deterministic and take about 1 h on
   a laptop. With `RESULTS_DIR=results_main`, `SIM_METHS=SC_true,SC_tight`, `SIM_SEEDS=0,1,2,3,4,5`,
   `SIM_CASES=1,2,3,4,5,6,7,8,9,10`, run `matlab/sim_run.m` once with `SIM_CUBE=pavia120`, `SIM_TAG=ctrlP` and once
   with `SIM_CUBE=belt100`, `SIM_TAG=ctrlB` (`SC_true`: the same model with the true per-band noise levels;
   `SC_tight`: tolerance 1e-7, at most 300 iterations). The predicted stripe false-alarm rate (Proposition 2) is
   computed by `python fa_predict.py > results_main/fa_predict.txt`.
5. `python make_rsl.py` (with `RESULTS_DIR=results_main`) turns the CSV files into the LaTeX tables and the numbers
   quoted in the text.

### Results included

`experiments/results_main/` holds the outputs of the full run reported in the manuscript (Intel Core i7-14700KF,
32 GB RAM, MATLAB R2022a, Windows 11; see `env_main.txt`): per-cell scores of all methods on the 116 simulated cubes
(`sim_*.csv`), real-data metrics and classification accuracy, timings, full-size scenes, rank sensitivity, t-CTV
tuning, the Python checks of the stripe false-alarm rate and of the settling of the noise levels, and the figures.
The released FastHyMix code stops with an index error in 5 of the 6 runs of Beltsville case 7; these cells are
missing in its results.

## License

MIT (see `LICENSE`). `experiments/matlab/s3ttv_cpu.m` and `experiments/matlab/prox_s3ttv_paged.m` are derived from
S3TTV, which is also MIT-licensed; see `THIRD_PARTY_NOTICES.md`.

## Citation

W. Wang, "Self-calibrated subspace robust PCA for parameter-free removal of mixed noise and line stripes in
hyperspectral images," manuscript submitted for publication, 2026.
