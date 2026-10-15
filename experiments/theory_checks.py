"""Numerical checks backing the propositions in the manuscript."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import norm
from common import *
from rcctv_np import rcctv, soft
from runner import cube

# ---- P1: min_s lam|s| + beta/2 (r-s)^2 == Huber(tau=lam/beta) (scaled by beta) --------------------
rng = np.random.default_rng(0)
lam, beta = 0.094, 3.1
r = rng.normal(0, 0.3, 200000)
s = soft(r, lam / beta)
val = lam * np.abs(s) + beta / 2 * (r - s) ** 2
tau = lam / beta
hub = np.where(np.abs(r) <= tau, beta / 2 * r ** 2, lam * np.abs(r) - lam ** 2 / (2 * beta))
print('P1 max |closed-form - Huber| =', float(np.abs(val - hub).max()))
M = N = 64; p = 99
print('P1 tau_k/sigma_k with lam=6/sqrt(MN):', 6 / np.sqrt(M * N) * (np.sqrt(M * N) + np.sqrt(p)) / (3 * np.sqrt(2)),
      ' = sqrt2(1+sqrt(p/MN)) =', np.sqrt(2) * (1 + np.sqrt(p / (M * N))))

# ---- P3: stripe false-alarm rate on stripe-free Gaussian data -----------------------------------
print('\nP3 false alarm (stripe-free data): fraction of nonzero stripe amplitudes a_jk')
print('  c   sigma   observed   2*Phi(-c)')
for name in ('pavia64', 'belt64', 'wdc1'):
    clean = cube(name)
    for sg in (0.05, 0.10):
        Y = clean + sg * np.random.default_rng(7).standard_normal(clean.shape)
        for c in (2.0, 3.0, 4.0):
            _, info = rcctv(Y, T_mode='pen', G_mode='band', cT=c, sig_online=True, return_info=True)
            a = info['T'].reshape(clean.shape, order='F')[0]           # (N,p)
            print(f'  {c:.0f}  {sg:.2f}   {np.mean(a != 0):.4f}     {2 * norm.sf(c):.4f}   [{name}]')

# ---- P4: scale fixed point at termination --------------------------------------------------------
print('\nP4 relative change of sigma_k over the last 3 iterations (max over bands) / beta_k*sigma_k')
for name, c in (('pavia64', 5), ('belt64', 8), ('wdc1', 3)):
    clean = cube(name)
    Y = make_noise(clean, c, 3)
    _, info = rcctv(Y, T_mode='pen', G_mode='band', cT=3.0, sig_online=True, return_info=True)
    h = np.array(info['sig_hist'])
    rel = np.abs(h[-1] - h[-2]) / h[-1]
    rel3 = np.abs(h[-1] - h[-4]) / h[-1]
    print(f'  {name} c{c}: iters={info["iters"]}  max rel change last step={rel.max():.2e}  last 3 steps={rel3.max():.2e}  '
          f'median sigma_hat={np.median(h[-1]):.4f}')
