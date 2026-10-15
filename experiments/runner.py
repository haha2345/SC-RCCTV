"""Parallel experiment runner.  usage: python runner.py <set> <out.csv> <arm,arm,...> [seeds]
set: dev (WDC x2, simu_indian; tuning only)  |  locked (pavia64, belt64 seed0 official inputs)
     | fresh (Pavia/Belt 64-crop clean, my generator, seeds 1..S)"""
import sys, os, csv, time, json
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.io as sio
from concurrent.futures import ProcessPoolExecutor
from common import *
from rcctv_np import rcctv, score

B = ROOT + 'code/baselines/'


def cube(name):
    if name.startswith('wdc'):
        W = sio.loadmat(B + 'NGMeet/data/WDC.mat')['WDC']
        o = {'wdc1': (0, 0), 'wdc2': (100, 100)}[name]
        return W[o[0]:o[0] + 64, o[1]:o[1] + 64, :100].astype(float)
    if name == 'indian':
        return sio.loadmat(B + 'LRTDTV/simu_indian.mat')['simu_indian'][:64, :64, :100].astype(float)
    if name in ('pavia64', 'belt64'):
        return sio.loadmat(f'{LOCK}{name}/clean.mat')['clean']
    raise KeyError(name)


# arm name -> kwargs (oracle arms need per-case beta/sigma, filled in below)
ARMS = {
    'official': dict(),
    'oracle': dict(T_mode='ball2', G_mode='peng', _oracle=True),
    'oracleG': dict(T_mode='none', G_mode='peng', _oracle=True),
    'blind_c3': dict(T_mode='pen', G_mode='band', cT=3.0),
    'blind_c1': dict(T_mode='pen', G_mode='band', cT=1.0),
    'blind_c6': dict(T_mode='pen', G_mode='band', cT=6.0),
    'blind_c10': dict(T_mode='pen', G_mode='band', cT=10.0),
    'blind_c3_on': dict(T_mode='pen', G_mode='band', cT=3.0, sig_online=True),
    'oracle_band': dict(T_mode='ball2', G_mode='peng', _oracle=True, _band=True),   # true per-band sigma_k
    'on_c2': dict(T_mode='pen', G_mode='band', cT=2.0, sig_online=True),
    'on_c4': dict(T_mode='pen', G_mode='band', cT=4.0, sig_online=True),
    'on_noG': dict(T_mode='pen', G_mode='none', cT=3.0),                 # ablation: stripe only
    'on_noT': dict(T_mode='none', G_mode='band', sig_online=True),       # ablation: dense only
    'on_scalar': dict(T_mode='pen', G_mode='band', cT=3.0, sig_online=True, scalar_sigma=True),
    'vh': dict(T_mode='pen', G_mode='band', cT=3.0, sig_online=True, stripe_dirs='vh'),
    'on_rk6': dict(T_mode='pen', G_mode='band', cT=3.0, sig_online=True, rk=6),
    'on_rk12': dict(T_mode='pen', G_mode='band', cT=3.0, sig_online=True, rk=12),
}


def job(a):
    name, c, seed, arm = a
    clean = cube(name)
    if seed == 0 and name in ('pavia64', 'belt64'):
        Y, _ = load_locked(name, c)
    else:
        Y = make_noise(clean, c, seed)
    if os.environ.get('TRANSPOSE') == '1':       # horizontal (row-constant) stripes
        Y = Y.transpose(1, 0, 2); clean = clean.transpose(1, 0, 2)
    kw = dict(ARMS[arm])
    band = kw.pop('_band', False)
    if kw.pop('_oracle', False):
        kw.update(beta=oracle_beta(c, clean.size), sigma=oracle_sigma(c, seed, clean.shape[2], band))
    t = time.time()
    X, info = rcctv(Y, return_info=True, **kw)
    ps, ss = score(X, clean)
    return dict(cube=name, case=c, seed=seed, arm=arm, mpsnr=ps, mssim=ss, iters=info['iters'], sec=time.time() - t)


if __name__ == '__main__':
    sets = {'dev': (['wdc1', 'wdc2', 'indian'], [1, 2]),
            'locked': (['pavia64', 'belt64'], [0]),
            'fresh': (['pavia64', 'belt64'], None)}
    S, out, arms = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
    cubes, seeds = sets[S]
    if S == 'fresh':
        seeds = list(range(1, int(sys.argv[4]) + 1))
    cases = [int(x) for x in os.environ.get('CASES', '1,2,3,4,5,6,7,8').split(',')]
    tasks = [(n, c, s, a) for n in cubes for c in cases for s in seeds for a in arms]
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        res = list(ex.map(job, tasks, chunksize=1))
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
        w.writeheader(); w.writerows(res)
    print('wrote', out, len(res))
