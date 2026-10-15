import os
import numpy as np, scipy.io as sio
# project root = parent of new_method_2026-10-12/ (was the absolute laptop path; relative since 2026-10-14 so the
# checks also run on the main machine)
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')).replace('\\', '/') + '/'
LOCK = ROOT + 'data/locked_nhsi/'
# cases 1-8: official eight cases; 9-10: NEW band-varying (non-i.i.d.) Gaussian sigma_k ~ U(.02,.18)
SIGMA = {1: .05, 2: .05, 3: .10, 4: 0., 5: .05, 6: .10, 7: .05, 8: .10, 9: None, 10: None}
SP = {1: 0, 2: .05, 3: .05, 4: 0, 5: 0, 6: 0, 7: .05, 8: .05, 9: 0, 10: .05}
STR = {1: 0, 2: 0, 3: 0, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 0, 10: 1}
RHO = {4: .98, 5: .95, 6: .95, 7: .9025, 8: .9025, 10: .9025}   # official stripe_radius = RHO*numel*r*I/2


def oracle_beta(c, numel):
    return RHO[c] * numel * 0.05 * 0.5 / 2 if STR[c] else 0.0


def sigma_vec(seed, p):
    return np.random.default_rng(1000 + seed).uniform(0.02, 0.18, p)


def oracle_sigma(c, seed, p, band=False):
    """Known-sigma oracle.  band=False gives the scalar the original method must use."""
    if c >= 9:
        sv = sigma_vec(seed, p)
        return sv if band else float(np.sqrt(np.mean(sv ** 2)))
    return SIGMA[c]


def load_locked(sc, c):
    Y = sio.loadmat(f'{LOCK}{sc}/case{c}_in.mat')['Nhsi']
    clean = sio.loadmat(f'{LOCK}{sc}/clean.mat')['clean']
    return Y, clean


def make_noise(clean, c, seed, I=0.5, rate=0.05):
    """Same recipe as the official 8 cases, reverse-engineered from the locked inputs:
    Gaussian sigma; 5% salt&pepper (0/1); stripes: 5% of columns per band, U(-I,I)."""
    rng = np.random.default_rng(seed)
    M, N, p = clean.shape
    Y = clean.copy()
    if c >= 9:
        Y = Y + sigma_vec(seed, p)[None, None, :] * rng.standard_normal(Y.shape)
    elif SIGMA[c] > 0:
        Y = Y + SIGMA[c] * rng.standard_normal(Y.shape)
    if SP[c] > 0:
        m = rng.random(Y.shape) < SP[c]
        Y[m] = rng.integers(0, 2, size=int(m.sum())).astype(float)
    if STR[c]:
        for k in range(p):
            cols = rng.random(N) < rate
            Y[:, cols, k] += rng.uniform(-I, I, size=int(cols.sum()))[None, :]
    return Y
