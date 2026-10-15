"""Dev-set inputs (WDC x2, simulated Indian Pines; 64x64x100) for tuning baselines that need it (t-CTV lambda).
Same cubes and generator as the dev set used for SC-RCCTV's own constants (runner.py); seed 1, cases 1-10.
Writes results/inputs_dev/<cube>_c<case>_s1.mat with Nhsi, clean."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scipy.io as sio
from common import make_noise
from runner import cube

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'inputs_dev')
os.makedirs(out, exist_ok=True)
for name in ('wdc1', 'wdc2', 'indian'):
    clean = cube(name)
    for c in range(1, 11):
        Y = make_noise(clean, c, 1)
        sio.savemat(os.path.join(out, f'{name}_c{c}_s1.mat'), {'Nhsi': Y, 'clean': clean})
print('dev inputs written to', out)
