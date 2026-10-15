"""Tables + number macros for the Remote Sensing Letters version (paper_rsl/), straight from the CSVs.

Every number quoted in the text is a macro \\N{key} defined in paper_rsl/numbers.tex; a key whose data are not
there yet (main-machine runs: S3TTV, idle timing, large scale) is simply not defined and prints a red TBD.
Run from new_method_2026-10-12/:  python make_rsl.py
env RESULTS_DIR (default results) selects the result folder, e.g. RESULTS_DIR=results_main for the full re-run on the
main machine; it must contain the CSVs (sim_*, timing*, large, real_*, rank_sens, tune_tctv, sigma_stats,
real_picks) and the text outputs of the Python checks (check_fa.txt, check_fp.txt).
"""
import glob, os, re
import numpy as np, pandas as pd

R = os.environ.get('RESULTS_DIR', 'results').rstrip('/\\') + '/'
OUT = 'paper_rsl/'
os.makedirs(OUT + 'tables', exist_ok=True)
df = pd.concat([pd.read_csv(f) for f in glob.glob(R + 'sim_*.csv') if 'smoke' not in f])
df = df[df.meth != 'tCTV'].drop_duplicates(['cube', 'case', 'seed', 'meth'], keep='last')   # tCTV = default lambda, superseded by tCTVd
CUBES = ['pavia120', 'belt100']
CN = {'pavia120': 'Pavia University', 'belt100': 'Beltsville'}
NAME = {'RCCTV': 'RCCTV', 'LRTDTV': 'LRTDTV', 'E3DTV': 'E-3DTV', 'tCTVd': 't-CTV', 'OLD': 'Oracle', 'SC': 'SC-RCCTV',
        'FastHyMix': 'FastHyMix', 'L1HyMixDe': 'L1HyMixDe', 'S3TTV': 'S3TTV', 'Input': 'Noisy input'}
GRP = {'gs': [1, 2, 3], 'st': [4, 5, 6, 7, 8], 'hv': [9, 10], 'all': list(range(1, 11))}
TBD = r'\tbd'
K = {}                                            # macro key -> text


def sgn(x, nd=2):
    s = f'{x:+.{nd}f}'
    if float(s) == 0:
        return f'{0:.{nd}f}'
    return s.replace('-', '$-$')


def fsec(x):
    return f'{x:.1f}' if x < 10 else f'{x:.0f}'


def bold_best(vals, higher=True, fmt='{:.2f}'):
    v = np.array(vals, float)
    if np.all(np.isnan(v)):
        return [TBD] * len(v)
    best = np.nanmax(v) if higher else np.nanmin(v)
    return [(TBD if np.isnan(x) else (r'\textbf{' + fmt.format(x) + '}' if x == best else fmt.format(x))) for x in v]


def save(name, lines):
    open(OUT + 'tables/' + name, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')


# ------------------------------------------------------------------ timing (idle machine; main-machine files only)
TM = None
if os.path.exists(R + 'timing.csv'):
    TM = pd.read_csv(R + 'timing.csv')[['ds', 'meth', 'sec']]
    TM = TM[TM.meth != 'tCTV']
    if os.path.exists(R + 'timing_new.csv'):
        TM = pd.concat([TM, pd.read_csv(R + 'timing_new.csv')[['ds', 'meth', 'sec']]])


def tsec(ds, m):
    if TM is None:
        return np.nan
    s = TM[(TM.ds == ds) & (TM.meth == m)].sec
    return float(s.iloc[0]) if len(s) else np.nan


_s3 = [tsec(d, 'S3TTV') for d in ('pavia120_c8', 'belt100_c8')]
if all(np.isfinite(_s3)):
    K['t.S3TTV.min'] = f'{np.mean(_s3) / 60:.0f}'                  # "about N min per cube"
S3MIN = K.get('t.S3TTV.min', TBD)

# ------------------------------------------------------------------ crashes of the released FastHyMix code
_sc_n = df[df.meth == 'SC'].groupby(['cube', 'case']).size()
_fh_n = df[df.meth == 'FastHyMix'].groupby(['cube', 'case']).size()
FHFAIL = [(c, k, int(n - _fh_n.get((c, k), 0)), int(n)) for (c, k), n in _sc_n.items()
          if 'FastHyMix' in set(df.meth) and _fh_n.get((c, k), 0) < n]
FH_TAIL = 'value of the remaining run' if FHFAIL and all(nt - nf == 1 for _, _, nf, nt in FHFAIL) else 'mean of the remaining runs'


def fh_note(tail=FH_TAIL):
    """footnote text on the FastHyMix crashes ('' if the released code never crashed)"""
    if not FHFAIL:
        return ''
    parts = [f'{nf} of the {nt} runs of {CN[c]} case~{k}' for c, k, nf, nt in FHFAIL]
    return 'FastHyMix stopped with an error inside the released code in ' + ' and '.join(parts) + '; ' + tail + '.'


K['fhfail.where'] = ' and '.join(f'{CN[c]} case~{k}' for c, k, _, _ in FHFAIL)

# ------------------------------------------------------------------ numbers from the auxiliary checks
def _read(name):
    return open(R + name, encoding='utf-8').read() if os.path.exists(R + name) else None


def sci(x):
    e = int(np.floor(np.log10(x)))
    return f'{x / 10 ** e:.1f}\\times10^{{{e}}}'


def fmt_bands(b):
    """[1, 61, 75, 88, 89, 90] -> '1, 61, 75 and 88--90'"""
    runs, start = [], None
    for i, x in enumerate(b):
        if start is None:
            start = x
        if i == len(b) - 1 or b[i + 1] != x + 1:
            runs.append(str(start) if start == x else f'{start}--{x}')
            start = None
    return runs[0] if len(runs) == 1 else ', '.join(runs[:-1]) + ' and ' + runs[-1]


_t = _read('check_fa.txt')                      # stripe false-alarm rate on stripe-free cubes (check_fa.py), in %
if _t:
    for c, nd in ((3, 2), (2, 1)):
        m = re.search(rf'c={c}: col range ([\d.]+)-([\d.]+), row range ([\d.]+)-([\d.]+)', _t)
        if m:
            v = [100 * float(x) for x in m.groups()]
            K[f'fa.c{c}'] = f'{min(v):.{nd}f}--{max(v):.{nd}f}'
_t = _read('fa_predict.txt')                     # predicted false-alarm rate (fa_predict.py), in %
if _t:
    for c, nd in ((3, 2), (2, 1)):
        m = re.search(rf'c={c}: predicted range ([\d.]+)-([\d.]+)%', _t)
        if m:
            lo, hi = f'{float(m.group(1)):.{nd}f}', f'{float(m.group(2)):.{nd}f}'
            K[f'fapred.c{c}'] = lo if lo == hi else f'{lo}--{hi}'
_t = _read('check_fp.txt')                      # settling of the noise levels at termination (check_fp.py)
if _t:
    m = re.search(r'runs:\s*(\d+).*?\(last step\):\s*([0-9.eE+-]+)\s+median over runs:\s*([0-9.eE+-]+)', _t)
    if m:
        K['fp.runs'], K['fp.max'], K['fp.med'] = m.group(1), sci(float(m.group(2))), sci(float(m.group(3)))
if os.path.exists(R + 'rank_sens.csv'):          # rank sensitivity on Indian Pines (rank_sens.m)
    _r = pd.read_csv(R + 'rank_sens.csv').set_index('rk')
    for rk in _r.index:
        K[f'rank.oa{rk}'] = f"{_r.loc[rk, 'OA']:.1f}"
_tf = glob.glob(R + 'tune_tctv*.csv')            # t-CTV weight multiplier chosen on the development cubes
if _tf:                                          # (one file per cube in the full re-run, tune_tctv_<cube>.csv)
    _u = pd.concat([pd.read_csv(f) for f in _tf]).drop_duplicates(['cube', 'case', 'mult'], keep='last')
    # paired comparison: only (cube, case) cells run with every multiplier
    _k = _u.groupby(['cube', 'case']).mult.nunique()
    _u = _u.set_index(['cube', 'case']).loc[_k[_k == _u.mult.nunique()].index].reset_index()
    K['tune.best'] = f"{_u.groupby('mult').mpsnr.mean().idxmax():g}"
if os.path.exists(R + 'sigma_stats.csv'):        # initial vs refined noise levels on Indian Pines (real_analyze.m)
    _s = pd.read_csv(R + 'sigma_stats.csv', dtype={'bands_up2': str}).iloc[0]
    K['sig.medratio'], K['sig.maxratio'] = f'{_s.med_ratio:.2f}', f'{_s.max_ratio:.1f}'
    K['sig.absdev'] = f'{100 * _s.abs_maxdev:.0f}'
    _b = [int(x) for x in str(_s.bands_up2).split()] if isinstance(_s.bands_up2, str) else []
    K['sig.bands2x'] = fmt_bands(_b) if _b else 'no band'
if os.path.exists(R + 'real_picks.csv'):         # bands shown in the real-data figure (real_analyze.m)
    _p = pd.read_csv(R + 'real_picks.csv').set_index('ds')
    K['pick.indian'], K['pick.urban'] = str(_p.loc['indian220', 'band']), str(_p.loc['urban200', 'band'])


for ds in ('pavia120_c8', 'belt100_c8', 'indian220'):
    for m in NAME:
        t = tsec(ds, m)
        if np.isfinite(t):
            K[f't.{m}.{ds}'] = fsec(t)
    sc = tsec(ds, 'SC')
    for m in NAME:
        t = tsec(ds, m)
        if np.isfinite(t) and np.isfinite(sc) and m != 'SC':
            r = t / sc
            K[f'tr.{m}.{ds}'] = f'{r:.1f}' if r < 10 else f'{r:.0f}'
# ratio ranges over the two native cubes
for m in NAME:
    rr = [tsec(d, m) / tsec(d, 'SC') for d in ('pavia120_c8', 'belt100_c8')]
    if all(np.isfinite(rr)):
        lo, hi = min(rr), max(rr)
        f = (lambda x: f'{x:.1f}') if hi < 10 else (lambda x: f'{x:.0f}')
        K[f'trr.{m}'] = f(lo) + '--' + f(hi) if f(lo) != f(hi) else f(lo)
        if lo >= 1000:                                   # "more than N times faster", N rounded down to hundreds
            K[f'trlo.{m}'] = f'{int(lo // 100 * 100):,}'          # Chicago style: 2,100

# ------------------------------------------------------------------ simulated: gains of SC over every reference
pv = df.pivot_table(index=['cube', 'case', 'seed'], columns='meth', values='mpsnr')
for ref in ['RCCTV', 'LRTDTV', 'E3DTV', 'FastHyMix', 'L1HyMixDe', 'OLD', 'OLDband', 'SCv', 'SC_scalar', 'SC_noT', 'SC_noG', 'SC_fix', 'SC_true', 'SC_tight', 'tCTVd', 'S3TTV']:
    if ref not in pv:
        continue
    d = (pv['SC'] - pv[ref]).dropna()
    if not len(d):
        continue
    K[f'min.{ref}'] = sgn(d.min())
    K[f'umin.{ref}'] = f'{abs(d.min()):.2f}'
    K[f'amax.{ref}'] = f'{d.abs().max():.2f}'           # largest |difference| in any cell
    for g, cs in GRP.items():
        dd = d[d.index.get_level_values('case').isin(cs)]
        if len(dd):
            K[f'd.{ref}.{g}'] = sgn(dd.mean())
            K[f'u.{ref}.{g}'] = f'{abs(dd.mean()):.2f}'          # unsigned, for "better/worse by"
            K[f'w.{ref}.{g}'] = f'{int((dd > 0).sum())} of {len(dd)}'
            K[f'n.{ref}.{g}'] = str(len(dd))
    for c in CUBES:
        dc = d.xs(c, level='cube')
        K[f'd.{ref}.{c}'] = sgn(dc.mean())
    for case, v in d.groupby(level='case').mean().items():
        K[f'd.{ref}.c{case}'] = sgn(v)
K['ncells'] = str(int(pv['SC'].notna().sum()))
K['fhfail'] = str(int(pv['SC'].notna().sum() - pv['FastHyMix'].notna().sum())) if 'FastHyMix' in pv else TBD
K['fhfail.sentence'] = (f"The released FastHyMix code stopped with an error in {K['fhfail']} runs ({K['fhfail.where']}); "
                        'these cells are excluded from its averages.') if FHFAIL else ''
# seed 0, cases 1-8 (the only cells with t-CTV and S3TTV): gains on the common cells
s0 = pv[(pv.index.get_level_values('seed') == 0) & (pv.index.get_level_values('case') <= 8)]
for ref in ['tCTVd', 'S3TTV', 'E3DTV', 'FastHyMix', 'L1HyMixDe', 'LRTDTV']:
    if ref not in s0:
        continue
    d = (s0['SC'] - s0[ref]).dropna()
    for g, cs in (('gs', [1, 2, 3]), ('st', [4, 5, 6, 7, 8]), ('all', list(range(1, 9)))):
        dd = d[d.index.get_level_values('case').isin(cs)]
        if len(dd):
            K[f's0.{ref}.{g}'] = sgn(dd.mean())
            K[f's0u.{ref}.{g}'] = f'{abs(dd.mean()):.2f}'
            K[f's0w.{ref}.{g}'] = f'{int((dd > 0).sum())} of {len(dd)}'
for m in pv.columns:
    for c in CUBES:
        v = s0.xs(c, level='cube')[m].dropna()
        if len(v) == 8:
            K[f'm0.{m}.{c}'] = f'{v.mean():.2f}'
# per cube, per case means
cm = df.groupby(['meth', 'cube', 'case']).mpsnr.mean()
for (m, c, case), v in cm.items():
    K[f'c.{m}.{c}.{case}'] = f'{v:.2f}'
# SC iterations and default-lambda t-CTV (for the text)
K['iters.SC'] = f"{df[df.meth == 'SC'].iters.mean():.1f}"
if 'S3TTV' in set(df.meth):
    it = df[df.meth == 'S3TTV'].iters
    K['s3.itmin'], K['s3.itmax'] = str(int(it.min())), str(int(it.max()))
old = pd.concat([pd.read_csv(f) for f in glob.glob(R + 'sim_tc[PB].csv')])
for c in CUBES:
    K[f'm0.tCTVdefault.{c}'] = f"{old[old.cube == c].mpsnr.mean():.2f}"

# ------------------------------------------------------------------ Table 1: MPSNR per case (two stacked panels)
# Only FastHyMix may have fewer cells than SC-RCCTV (released code crashed; footnote a). Any other method with
# missing cells is still running: its cell prints TBD instead of a mean over a partial set of realizations.
meths = ['RCCTV', 'LRTDTV', 'E3DTV', 'tCTVd', 'FastHyMix', 'L1HyMixDe', 'SC']
hdr = ['RCCTV', 'LRTDTV', 'E-3DTV', 't-CTV', 'FastHyMix', 'L1HyMixDe', 'Ours']
lines = [r'\begin{table}[!htbp]', r'\tbl{MPSNR (dB) on the simulated cubes, mean over six noise realizations (seed 0 = official inputs, cases 1--8; seeds 1--5, all cases). Best per row in bold.}',
         r'{\small\setlength{\tabcolsep}{4.2pt}', r'\begin{tabular}{c' + 'c' * len(meths) + '}', r'\toprule',
         'Case & ' + ' & '.join(hdr) + r'\\']
cnt = df[df.meth.isin(meths)].groupby(['cube', 'case', 'meth']).size()
for c in CUBES:
    sizes = {'pavia120': r'$120{\times}120{\times}99$', 'belt100': r'$100{\times}100{\times}128$'}
    lines += [r'\midrule', r'\multicolumn{%d}{l}{\emph{%s} (%s)}\\' % (len(meths) + 1, CN[c], sizes[c])]
    for case in range(1, 11):
        full = cnt.get((c, case, 'SC'), 0)
        v = [cm.get((m, c, case), np.nan) if (m == 'FastHyMix' or cnt.get((c, case, m), 0) == full) else np.nan for m in meths]
        b = bold_best(v)
        b = [(r'$\times$' if m == 'FastHyMix' and cnt.get((c, case, m), 0) == 0 and 'FastHyMix' in set(df.meth) else x)
             + (r'$^{\rm a}$' if m == 'FastHyMix' and 0 < cnt.get((c, case, m), 0) < full else '') for x, m in zip(b, meths)]
        lines.append(' & '.join([str(case)] + b) + r'\\')
lines += [r'\bottomrule', r'\end{tabular}}'] + ([r'\tabnote{$^{\rm a}$' + fh_note() + '}'] if FHFAIL else []) + [
          r'\label{tab:mpsnr}', r'\end{table}']
save('tab_mpsnr.tex', lines)

# ------------------------------------------------------------------ Table 2: all metrics averaged + time
allm = ['RCCTV', 'LRTDTV', 'E3DTV', 'tCTVd', 'FastHyMix', 'L1HyMixDe', 'SC']
lines = [r'\begin{table}[!htbp]', r'\tbl{Mean over all cases and realizations of the simulated cubes. Time: one run on case~8 (seed~0) on an otherwise idle machine. Below the inner rule: S3TTV, run on seed~0, cases 1--8 only (about ' + S3MIN + r'\,min per cube), compared with SC-RCCTV on the same 8 inputs (per-case values in Table~S1 of the supplementary material).}',
         r'{\small\setlength{\tabcolsep}{5pt}', r'\begin{tabular}{llccccc}', r'\toprule',
         r'Cube & Method & MPSNR$\uparrow$ & MSSIM$\uparrow$ & SAM$\downarrow$ & ERGAS$\downarrow$ & Time (s)\\']
s0df = df[(df.seed == 0) & (df.case <= 8)]
for c in CUBES:
    sub = df[(df.cube == c) & df.meth.isin(allm)]
    full = sub[sub.meth == 'SC'][['case', 'seed']]
    g = sub.groupby('meth')[['mpsnr', 'mssim', 'sam', 'ergas']].mean()
    lines.append(r'\midrule')
    complete = [m for m in allm if m in g.index and (m == 'FastHyMix' or len(sub[sub.meth == m]) == len(full))]
    gb = g.loc[complete]                                    # best values only among methods with all cells
    for i, m in enumerate(allm):
        t = tsec(c + '_c8', m)
        ts = TBD if np.isnan(t) else fsec(t)
        lead = CN[c].split()[0] if i == 0 else ''
        if m not in complete:                               # still running: no mean over a partial set of cells
            lines.append(f"{lead} & {NAME[m]} & {' & '.join([TBD] * 4)} & {ts}" + r'\\')
            continue
        r = g.loc[m]
        f = lambda x, b, fm: (r'\textbf{' + fm.format(x) + '}') if x == b else fm.format(x)
        dag = r'$^{\rm a}$' if len(sub[sub.meth == m]) < len(full) else ''
        lines.append(f"{lead} & {NAME[m]}{dag} & {f(r.mpsnr, gb.mpsnr.max(), '{:.2f}')} & {f(r.mssim, gb.mssim.max(), '{:.4f}')} & "
                     f"{f(r.sam, gb.sam.min(), '{:.2f}')} & {f(r.ergas, gb.ergas.min(), '{:.2f}')} & {ts}" + r'\\')
        K[f'avg.{m}.{c}'] = f'{r.mpsnr:.2f}'
        K[f'avgssim.{m}.{c}'] = f'{r.mssim:.3f}'
        K[f'avgsam.{m}.{c}'] = f'{r.sam:.2f}'
    # seed 0, cases 1-8: S3TTV against SC-RCCTV on the same inputs
    lines.append(r'\cmidrule{2-7}')
    sub0 = s0df[s0df.cube == c]
    g0 = sub0.groupby('meth')[['mpsnr', 'mssim', 'sam', 'ergas']].mean()
    n0 = sub0.groupby('meth').size()
    m0 = [m for m in ('S3TTV', 'SC') if m in g0.index and n0[m] == 8]
    for m in ('S3TTV', 'SC'):
        t = tsec(c + '_c8', m)
        ts = TBD if np.isnan(t) else fsec(t)
        if m in m0:
            r = g0.loc[m]
            f = lambda x, col, fm, hi: (r'\textbf{' + fm.format(x) + '}') if x == (g0.loc[m0, col].max() if hi else g0.loc[m0, col].min()) else fm.format(x)
            vals = f"{f(r.mpsnr, 'mpsnr', '{:.2f}', 1)} & {f(r.mssim, 'mssim', '{:.4f}', 1)} & {f(r.sam, 'sam', '{:.2f}', 0)} & {f(r.ergas, 'ergas', '{:.2f}', 0)}"
        else:
            vals = ' & '.join([TBD] * 4)
        lines.append(f" & {NAME[m]}$^{{\\rm b}}$ & {vals} & {ts}" + r'\\')
lines += [r'\bottomrule', r'\end{tabular}}', r'\tabnote{' + (r'$^{\rm a}$Mean over the cells in which FastHyMix ran (Table~\ref{tab:mpsnr}). ' if FHFAIL else '') + r'$^{\rm b}$Seed~0, cases 1--8 only.}', r'\label{tab:avg}', r'\end{table}']
save('tab_avg.tex', lines)

# ------------------------------------------------------------------ Table 3: seed 0, cases 1-8 with t-CTV and S3TTV
m3 = ['E3DTV', 'FastHyMix', 'tCTVd', 'S3TTV', 'SC']
h3 = ['E-3DTV', 'FastHyMix', 't-CTV', 'S3TTV', 'Ours']
lines = [r'\begin{table}[!htbp]', r'\tbl{MPSNR (dB) on the official inputs (seed~0, cases 1--8), the only cells on which S3TTV was run (about ' + S3MIN + r'\,min per cube). t-CTV: weight tuned on development cubes (Section~5.1 of the letter).}',
         r'{\small\setlength{\tabcolsep}{3.4pt}', r'\begin{tabular}{c' + 'c' * len(m3) + '|' + 'c' * len(m3) + '}', r'\toprule',
         r' & \multicolumn{%d}{c|}{Pavia University} & \multicolumn{%d}{c}{Beltsville}\\' % (len(m3), len(m3)),
         'Case & ' + ' & '.join(h3) + ' & ' + ' & '.join(h3) + r'\\', r'\midrule']
ran = set(df.meth)                                   # a method with rows but a missing cell = failed run, not pending
failed = False


def mark_failed(cells, v):
    global failed
    out = []
    for x, m, val in zip(cells, m3, v):
        if np.isnan(val) and m in ran:
            out.append(r'$\times$'); failed = True
        else:
            out.append(x)
    return out


for case in range(1, 9):
    row = [str(case)]
    for c in CUBES:
        v = [pv.loc[(c, case, 0), m] if (m in pv and (c, case, 0) in pv.index) else np.nan for m in m3]
        row += mark_failed(bold_best(v), v)
    lines.append(' & '.join(row) + r'\\')
lines.append(r'\midrule')
row = ['Mean']
for c in CUBES:
    v = [s0.xs(c, level='cube')[m].mean() if m in s0 and s0.xs(c, level='cube')[m].notna().sum() == 8 else np.nan for m in m3]
    row += [('--' if (np.isnan(val) and m in ran) else x) for x, m, val in zip(bold_best(v), m3, v)]
lines.append(' & '.join(row) + r'\\')
row = ['Time (s)']
for c in CUBES:
    row += [TBD if np.isnan(tsec(c + '_c8', m)) else fsec(tsec(c + '_c8', m)) for m in m3]
lines.append(' & '.join(row) + r'\\')
lines += [r'\bottomrule', r'\end{tabular}}']
if failed:
    lines.append(r'\tabnote{$\times$: the released code stopped with an error; no mean is given for that method.}')
lines += [r'\label{tab:seed0}', r'\end{table}']
save('tab_seed0.tex', lines)

# ------------------------------------------------------------------ Table 4: ablation
fam = {'Gauss.+imp. (1--3)': [1, 2, 3], 'Stripes (4--8)': [4, 5, 6, 7, 8], 'Band-var. (9--10)': [9, 10]}
variants = [('SC', 'SC-RCCTV (full)'), ('SC_noT', 'without stripe term $T$'), ('SC_noG', 'without dense term $G$'),
            ('SC_scalar', 'scalar $\\sigma$ instead of $\\sigma_k$'),
            ('SC_fix', '$\\hat\\sigma_k^{(0)}$ kept fixed (two-stage)'), ('SCv', 'column stripes only'),
            ('SC_true', 'true $\\sigma_k$ (all else unchanged)'),
            ('OLD', 'RCCTV-K, true $\\sigma$ and radius'), ('OLDband', 'RCCTV-K, true $\\sigma_k$ and radius')]
known = (' Below the rule: SC-RCCTV with the true $\\sigma_k$ in place of the self-calibrated ones, all else unchanged;'
         ' and RCCTV-K, a different noise-aware model: RCCTV with a dense term weighted by the true noise level (scalar $\\sigma$,'
         ' the root mean square of $\\sigma_k$ in cases 9--10, or the band-wise $\\sigma_k$) and an $\\ell_1$-ball stripe set with the true radius.')
for fname, keep, label in (('tab_abl.tex', ['SC', 'SC_noT', 'SC_noG', 'SC_scalar', 'SC_fix', 'SC_true', 'OLD', 'OLDband'], 'tab:abl'),
                           ('tab_abl_full.tex', [k for k, _ in variants], 'tab:ablfull')):
    lines = [r'\begin{table}[!htbp]', r'\tbl{Ablation: mean MPSNR gain (dB) over RCCTV by noise family (both cubes, all realizations).' + known + '}',
             r'{\small\setlength{\tabcolsep}{6pt}', r'\begin{tabular}{lcccc}', r'\toprule',
             r'Variant & ' + ' & '.join(fam) + r' & All\\', r'\midrule']
    for key, lab in variants:
        if key not in keep:
            continue
        if key == 'SC_true':
            lines.append(r'\midrule')
        if key not in pv:
            lines.append(lab + ' & ' + ' & '.join([TBD] * (len(fam) + 1)) + r'\\')
            continue
        d = (pv[key] - pv['RCCTV']).dropna()
        cells = [sgn(d[d.index.get_level_values('case').isin(cs)].mean()) for cs in fam.values()] + [sgn(d.mean())]
        lines.append(lab + ' & ' + ' & '.join(cells) + r'\\')
    lines += [r'\bottomrule', r'\end{tabular}}', r'\label{' + label + '}', r'\end{table}']
    save(fname, lines)
# online vs fixed noise levels: other indices and the largest gain (for the text)
if 'SC_fix' in pv:
    for col, nd in (('mssim', 3), ('sam', 2)):                 # per cube, SC-RCCTV minus the fixed-sigma variant
        pc = df.pivot_table(index=['cube', 'case', 'seed'], columns='meth', values=col)
        dd = (pc['SC'] - pc['SC_fix']).dropna()
        for c in CUBES:
            K[f'dfix.{col}.{c}'] = sgn(dd.xs(c, level='cube').mean(), nd)
            K[f'ufix.{col}.{c}'] = f'{abs(dd.xs(c, level="cube").mean()):.{nd}f}'
    K['ufix.mssim.max'] = f"{max(float(K[f'ufix.mssim.{c}']) for c in CUBES):.3f}"
    d = (pv['SC'] - pv['SC_fix']).dropna()
    K['fix.max'] = f'{d.max():.2f}'
    K['fix.c4'] = f"{d[d.index.get_level_values('case') == 4].mean():.2f}"
    _o = d[d.index.get_level_values('case') != 4]
    K['fix.other'] = f'{_o.mean():.2f}'                                    # cases other than 4
    K['fix.share'] = f"{100 * d[d.index.get_level_values('case') == 4].sum() / d.sum():.0f}"   # % of the total gain from case 4
# stopping-rule check: tolerance 1e-7 and up to 300 iterations instead of 1e-5 / 100
if 'SC_tight' in pv:
    _d = (pv['SC_tight'] - pv['SC']).dropna().abs()
    K['tight.max'], K['tight.mean'], K['tight.n'] = f'{_d.max():.2f}', f'{_d.mean():.2f}', str(len(_d))
    K['iters.SC_tight'] = f"{df[df.meth == 'SC_tight'].iters.mean():.0f}"

# ------------------------------------------------------------------ Table 5: full-size scenes
LM = ['RCCTV', 'LRTDTV', 'E3DTV', 'tCTVd', 'FastHyMix', 'L1HyMixDe', 'SC']
L = pd.read_csv(R + 'large.csv') if os.path.exists(R + 'large.csv') else None
lines = [r'\begin{table}[!htbp]', r'\tbl{Full-size scenes, one run per method on an otherwise idle machine: Pavia University $610{\times}340{\times}103$ (cases 8 and 10, one realization) and Urban $307{\times}307{\times}210$ (real; stripe index averaged over the 20 bands with the highest stripe index in the input, $\times10^{-3}$).}',
         r'{\small\setlength{\tabcolsep}{5pt}', r'\begin{tabular}{lccccc}', r'\toprule',
         r' & \multicolumn{3}{c}{Pavia University (full)} & \multicolumn{2}{c}{Urban (full)}\\',
         r'Method & MPSNR c8 & MPSNR c10 & Time (s) & Stripe$\downarrow$ & Time (s)\\', r'\midrule']


def lget(job, m, col):
    if L is None:
        return np.nan
    s = L[(L.job == job) & (L.meth == m)][col]
    return float(s.iloc[0]) if len(s) else np.nan


p8 = [lget('paviaU_c8', m, 'mpsnr') for m in LM]
p10 = [lget('paviaU_c10', m, 'mpsnr') for m in LM]
ust = [lget('urban307', m, 'stripe_top20') * 1e3 for m in LM]
b8, b10, bst = bold_best(p8), bold_best(p10), bold_best(ust, higher=False, fmt='{:.1f}')
for i, m in enumerate(LM):
    tp = np.nanmean([lget('paviaU_c8', m, 'sec'), lget('paviaU_c10', m, 'sec')]) if L is not None else np.nan
    tu = lget('urban307', m, 'sec')
    lines.append(f"{NAME[m]} & {b8[i]} & {b10[i]} & {TBD if np.isnan(tp) else fsec(tp)} & {bst[i]} & {TBD if np.isnan(tu) else fsec(tu)}" + r'\\')
    for k, v in ((f'L.{m}.c8', p8[i]), (f'L.{m}.c10', p10[i]), (f'L.{m}.urbanst', ust[i])):
        if np.isfinite(v):
            K[k] = f'{v:.2f}' if 'urban' not in k else f'{v:.1f}'
    if np.isfinite(tp):
        K[f'L.{m}.tpavia'] = fsec(tp)
    if np.isfinite(tu):
        K[f'L.{m}.turban'] = fsec(tu)
lines += [r'\bottomrule', r'\end{tabular}}',
          r'\tabnote{S3TTV is not included: its memory use exceeds 16\,GB at this size.}',
          r'\label{tab:large}', r'\end{table}']
save('tab_large.tex', lines)
if L is not None:
    sc = [lget(j, 'SC', 'sec') for j in ('paviaU_c8', 'paviaU_c10', 'urban307')]
    for m in LM:
        if m == 'SC':
            continue
        rr = [lget(j, m, 'sec') / s for j, s in zip(('paviaU_c8', 'paviaU_c10', 'urban307'), sc)]
        if all(np.isfinite(rr)):
            K[f'Ltr.{m}'] = f'{np.mean(rr):.1f}' if np.mean(rr) < 10 else f'{np.mean(rr):.0f}'

# ------------------------------------------------------------------ Table 6: real data
rm = pd.read_csv(R + 'real_metrics.csv')
cl = pd.read_csv(R + 'real_indian_cls.csv').set_index('meth')
order = ['Input', 'RCCTV', 'LRTDTV', 'E3DTV', 'tCTVd', 'FastHyMix', 'L1HyMixDe', 'SC']
lines = [r'\begin{table}[!htbp]', r'\tbl{Real data. Indian Pines: overall accuracy (OA, mean $\pm$ std over 5 splits) and $\kappa$ of an RBF-SVM with 10\% training samples per class. Urban ($200{\times}200$ crop): stripe index averaged over the 20 noisiest bands of the input, and largest change of a band mean, both $\times10^{-3}$. Time: Indian Pines, idle machine; t-CTV was run on Indian Pines only.}',
         r'{\small\setlength{\tabcolsep}{5pt}', r'\begin{tabular}{lccccc}', r'\toprule',
         r'Method & OA (\%)$\uparrow$ & $\kappa\uparrow$ & Stripe$\downarrow$ & Shift$\downarrow$ & Time (s)\\', r'\midrule']
for k in order:
    u = rm[(rm.ds == 'urban200') & (rm.meth == k)]
    oa = f"{cl.loc[k, 'OA']:.1f}\\,$\\pm$\\,{cl.loc[k, 'OA_std']:.1f}" if k in cl.index else (TBD if k == 'tCTVd' else '--')
    kp = f"{cl.loc[k, 'kappa']:.3f}" if k in cl.index else (TBD if k == 'tCTVd' else '--')
    si = f"{1e3 * u.stripeIdx_top20.iloc[0]:.1f}" if len(u) else '--'
    sh = ('0' if k == 'Input' else f"{1e3 * u.shift_max.iloc[0]:.1f}") if len(u) else '--'
    t = tsec('indian220', k)
    sec = '--' if k == 'Input' else (TBD if np.isnan(t) else fsec(t))
    lines.append(f'{NAME[k]} & {oa} & {kp} & {si} & {sh} & {sec}' + r'\\')
    if k in cl.index:
        K[f'oa.{k}'] = f"{cl.loc[k, 'OA']:.1f}"
    if len(u):
        K[f'si.{k}'] = f"{1e3 * u.stripeIdx_top20.iloc[0]:.1f}"
        K[f'sh.{k}'] = f"{1e3 * u.shift_max.iloc[0]:.1f}"
lines += [r'\bottomrule', r'\end{tabular}}', r'\label{tab:real}', r'\end{table}']
save('tab_real.tex', lines)
for k in ('SC_noT', 'SC_noG', 'SC_scalar', 'SC_fix', 'tCTV'):
    if k in cl.index:
        K[f'oa.{k}'] = f"{cl.loc[k, 'OA']:.1f}"
for k in ('SC_noT', 'SC_fix'):
    u = rm[(rm.ds == 'urban200') & (rm.meth == k)]
    if len(u):
        K[f'si.{k}'] = f"{1e3 * u.stripeIdx_top20.iloc[0]:.1f}"

# ------------------------------------------------------------------ numbers.tex
with open(OUT + 'numbers.tex', 'w', encoding='utf-8') as fh:
    fh.write('% generated by make_rsl.py -- do not edit\n')
    for k in sorted(K):
        fh.write('\\expandafter\\def\\csname n@%s\\endcsname{%s}\n' % (k, K[k]))
missing_hint = [k for k in ('t.SC.pavia120_c8', 's0.S3TTV.all', 'L.SC.c8', 'oa.tCTVd') if k not in K]
print(f'{len(K)} macros written; not yet available: {missing_hint}')
