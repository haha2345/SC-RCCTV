"""Predicted stripe false-alarm rate of SC-RCCTV (Proposition 2 of the manuscript) for the cubes of check_fa.py.

At a stationary point with fixed weights, a stripe-free column receives a zero amplitude iff the mean of the clipped
Huber scores psi(r) = clip(beta*r, -lambda, lambda) of its residuals lies within +-lambda_v. For i.i.d. N(0, sigma^2)
residuals the central limit theorem gives the activation probability 2*Phi(-c / sqrt(kappa(t))), with t = tau/sigma
and kappa(t) = E[clip(Z, -t, t)^2] < 1, i.e. below the nominal 2*Phi(-c). A Monte Carlo of the exact column problem
checks the approximation.  Usage (from this folder):  python fa_predict.py > results_main/fa_predict.txt
"""
import numpy as np
from scipy.stats import norm


def kappa(t):
    return (2 * norm.cdf(t) - 1) - 2 * t * norm.pdf(t) + 2 * t * t * norm.sf(t)


def rate(c, t):
    return 2 * norm.sf(c / np.sqrt(kappa(t)))


cubes = {'pavia64': (64, 64, 99), 'belt64': (64, 64, 128), 'wdc1': (64, 64, 100)}   # sizes of the check_fa.py cubes
rng = np.random.default_rng(0)
print('cube     t=tau/sigma  kappa   predicted c=3   predicted c=2   Monte Carlo c=3   Monte Carlo c=2')
res = {3: [], 2: []}
for name, (M, N, p) in cubes.items():
    t = np.sqrt(2) * (1 + np.sqrt(p / (M * N)))
    mc = {}
    for c in (3, 2):
        r = np.clip(rng.standard_normal((200000, M)), -t, t)          # beta*sigma = 1 units
        mc[c] = float(np.mean(np.abs(r.mean(1)) > c / np.sqrt(M)))
        res[c].append(rate(c, t))
    print(f'{name:8s} {t:9.3f}  {kappa(t):6.3f}   {100 * rate(3, t):10.3f}%   {100 * rate(2, t):11.2f}%   '
          f'{100 * mc[3]:12.3f}%   {100 * mc[2]:13.2f}%')
for c in (3, 2):
    print(f'c={c}: predicted range {100 * min(res[c]):.3f}-{100 * max(res[c]):.3f}%  (nominal 2Phi(-c) {100 * 2 * norm.sf(c):.3f}%)')
