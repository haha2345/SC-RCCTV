# Third-party code (not included)

`experiments/matlab/setup_paths.m` expects these folders here. Obtain them from the authors:

| folder | method | source |
|---|---|---|
| `official_rcctv/` (with `TV_operator/`, `PROPACK/`) | RCCTV (`RCCTV_RPCA.m`) | https://github.com/mendy-2013/RCCTV |
| `official_tctv/` (with `lib/`, `TCTV/`) | t-CTV (`TCTV_TRPCA.m`) | https://github.com/wanghailin97 (Guaranteed Tensor Recovery Fused Low-rankness and Smoothness) |
| `baselines/LRTDTV/` (with `prox_operators/`, `mylib/`, `simu_indian.mat`) | LRTDTV | https://github.com/zhaoxile/Hyperspectral-Image-Restoration-via-Total-Variation-Regularized-Low-rank-Tensor-Decomposition |
| `baselines/NGMeet/data/WDC.mat` | Washington DC Mall cube (development data only) | https://github.com/quanmingyao/NGMeet |

After placing t-CTV, create the variant used in the experiments:
`python tools/make_tctv_trpca_lam.py code/official_tctv/TCTV_TRPCA.m`.
