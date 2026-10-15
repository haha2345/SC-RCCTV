"""False-alarm rate of the final model (stripe_dirs='vh') on stripe-free Gaussian cubes."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import norm
from common import *
from rcctv_np import rcctv
from runner import cube

print('  c  sigma  col-amp(a)  row-amp(b)   nominal 2Phi(-c)')
res = {2.0: [], 3.0: []}
for name in ('pavia64', 'belt64', 'wdc1'):
    clean = cube(name)
    M, N, p = clean.shape
    for sg in (0.05, 0.10):
        Y = clean + sg * np.random.default_rng(7).standard_normal(clean.shape)
        for c in (2.0, 3.0):
            _, info = rcctv(Y, T_mode='pen', G_mode='band', cT=c, sig_online=True, stripe_dirs='vh', return_info=True)
            fa_a = float(np.mean(info['a'] != 0))      # column-stripe amplitudes a_jk
            fa_b = float(np.mean(info['b'] != 0))      # row-stripe amplitudes b_ik
            res[c].append((fa_a, fa_b))
            print(f'  {c:.0f}  {sg:.2f}    {fa_a:.4f}      {fa_b:.4f}       {2 * norm.sf(c):.4f}   [{name}]')
for c in res:
    r = np.array(res[c]); print(f'c={c:.0f}: col range {r[:, 0].min():.4f}-{r[:, 0].max():.4f}, row range {r[:, 1].min():.4f}-{r[:, 1].max():.4f}')
