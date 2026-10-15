"""Write native-size simulation inputs for the MATLAB driver.
cubes: pavia120 (120x120x99), belt100 (100x100x128).  seed 0 = official locked inputs (cases 1-8);
seeds 1..S = my generator (cases 1-10; 9/10 = band-varying sigma_k).  Each file: Nhsi, clean-free meta
(sigma_true 1xp, beta_true scalar, case, seed).  Clean cube is read from data/locked_nhsi/<cube>/clean.mat."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.io as sio
from common import *

S = int(sys.argv[1]) if len(sys.argv) > 1 else 5
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'inputs')
os.makedirs(out, exist_ok=True)
for cube in ('pavia120', 'belt100'):
    clean = sio.loadmat(f'{LOCK}{cube}/clean.mat')['clean']
    p = clean.shape[2]
    for c in range(1, 11):
        for s in range(0, S + 1):
            if s == 0 and c > 8:
                continue
            Y = load_locked(cube, c)[0] if s == 0 else make_noise(clean, c, s)
            sg = oracle_sigma(c, s, p, band=True)
            sg = np.full(p, sg, float) if np.isscalar(sg) else np.asarray(sg, float)
            sio.savemat(f'{out}/{cube}_c{c}_s{s}.mat',
                        dict(Nhsi=Y, sigma_true=sg[None, :], beta_true=oracle_beta(c, clean.size), case=c, seed=s))
    print(cube, 'done')
