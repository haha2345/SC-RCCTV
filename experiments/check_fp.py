"""Scale fixed-point check at termination over many runs (final model: stripe_dirs='vh')."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import *
from rcctv_np import rcctv
from runner import cube

rows = []
for name in ('pavia64', 'belt64'):
    clean = cube(name)
    for c in range(1, 11):
        for s in (0, 1):
            if s == 0 and c > 8:
                continue
            Y = load_locked(name, c)[0] if s == 0 else make_noise(clean, c, s)
            _, info = rcctv(Y, T_mode='pen', G_mode='band', sig_online=True, stripe_dirs='vh', return_info=True)
            h = np.array(info['sig_hist'])
            ok = h[-1] > 1e-4
            rel = np.abs(h[-1] - h[-2])[ok] / h[-1][ok]
            rows.append((name, c, s, info['iters'], float(rel.max()), float(np.median(rel))))
r = np.array([x[4] for x in rows])
print('runs:', len(rows), ' max over runs of max_k rel change (last step): %.2e   median over runs: %.2e' % (r.max(), np.median(r)))
print('iterations: min %d max %d mean %.1f' % (min(x[3] for x in rows), max(x[3] for x in rows), np.mean([x[3] for x in rows])))
worst = max(rows, key=lambda x: x[4]); print('worst run:', worst)
