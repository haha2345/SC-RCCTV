"""Large-scale simulated inputs: full Pavia University (610x340x103), global min-max to [0,1],
cases 8 (sigma=0.1 + 5% impulses + stripes) and 10 (band-varying sigma_k + impulses + stripes), seed 1,
same generator as the native-size experiment (common.make_noise).  -> results/inputs_large/"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.io as sio
from common import make_noise, oracle_sigma

here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, '..', 'data', 'large', 'PaviaU.mat')
out = os.path.join(here, 'results', 'inputs_large'); os.makedirs(out, exist_ok=True)
X = sio.loadmat(src)['paviaU'].astype(float)
X = (X - X.min()) / (X.max() - X.min())
for c in (8, 10):
    Y = make_noise(X, c, 1)
    sig = oracle_sigma(c, 1, X.shape[2], band=True)
    sio.savemat(os.path.join(out, f'paviaU_c{c}_s1.mat'), {'Nhsi': Y, 'clean': X, 'sigma_true': np.atleast_1d(sig)}, do_compression=False)
    print(c, Y.shape, float(Y.min()), float(Y.max()))
