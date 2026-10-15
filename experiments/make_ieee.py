"""IEEE-format tables for the GRSL letter (paper_grsl/) and the JSTARS article (paper_jstars/).

Runs make_rsl.py first, in this namespace, so every number in the three versions comes from the same data and
the same macro definitions; numbers.tex is copied from paper_rsl/. Cells of a method whose runs are not all in
results/ yet print TBD (only FastHyMix may legitimately miss cells: its released code crashed, footnote dagger).
Run from new_method_2026-10-12/:  python make_ieee.py
"""
import os
import shutil

import numpy as np

_ns = {'__name__': 'make_rsl'}
exec(compile(open('make_rsl.py', encoding='utf-8').read(), 'make_rsl.py', 'exec'), _ns)
df, pv, cm, K, TBD = _ns['df'], _ns['pv'], _ns['cm'], _ns['K'], _ns['TBD']
FHFAIL, fh_note, FH_TAIL, S3MIN = _ns['FHFAIL'], _ns['fh_note'], _ns['FH_TAIL'], _ns['S3MIN']
tsec, fsec, sgn, bold_best = _ns['tsec'], _ns['fsec'], _ns['sgn'], _ns['bold_best']
CUBES, CN, NAME, rm, cl, lget = _ns['CUBES'], _ns['CN'], _ns['NAME'], _ns['rm'], _ns['cl'], _ns['lget']

GRSL, JSTARS = 'paper_grsl/', 'paper_jstars/'
BOTH = [GRSL, JSTARS]
SHORT = {'RCCTV': 'RCCTV', 'LRTDTV': 'LRTD', 'E3DTV': 'E3DTV', 'tCTVd': 't-CTV', 'FastHyMix': 'FHMix',
         'L1HyMixDe': 'L1HMD', 'S3TTV': 'S3TTV', 'SC': 'Ours'}
LEGEND = 'LRTD: LRTDTV; FHMix: FastHyMix; L1HMD: L1HyMixDe.'
SIZE = {'pavia120': r'$120{\times}120{\times}99$', 'belt100': r'$100{\times}100{\times}128$'}
FAM = {'G+I': [1, 2, 3], 'Str.': [4, 5, 6, 7, 8], 'B.-v.': [9, 10]}
M7 = ['RCCTV', 'LRTDTV', 'E3DTV', 'tCTVd', 'FastHyMix', 'L1HyMixDe', 'SC']
DAG = r'$^\dagger$'
# footnotes on the crashes of the released FastHyMix code, from the data ('' when it never crashed)
DAGNOTE_CASE = (DAG + fh_note(FH_TAIL)) if FHFAIL else ''
DAGNOTE = (DAG + fh_note('means over the remaining cells')) if FHFAIL else ''
FAILCASES = {c: {k for cc, k, _, _ in FHFAIL if cc == c} for c in CUBES}


def write(name, lines, dirs=BOTH):
    for d in dirs:
        os.makedirs(d + 'tables', exist_ok=True)
        open(d + 'tables/' + name, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')


def note(text):
    return [r'\par\vspace{2pt}', r'{\raggedright\scriptsize ' + text + r'\par}']


cnt = df[df.meth.isin(M7 + ['S3TTV'])].groupby(['cube', 'case', 'meth']).size()
ncell = df[df.meth == 'SC'].groupby('cube').size()


def complete(m, c):
    """all cells of cube c present (FastHyMix: at least one, its missing cells are the documented crashes)"""
    n = int(df[(df.cube == c) & (df.meth == m)].shape[0])
    return n > 0 if m == 'FastHyMix' else n == int(ncell[c])


# ------------------------------------------------------------------ per-case MPSNR, stacked cubes (one column)
lines = [r'\begin{table}[!t]',
         r'\caption{MPSNR (dB) on the simulated cubes, mean over six noise realizations (seed~0: official inputs, cases 1--8; seeds 1--5: all cases). Best per row in bold. ' + LEGEND + '}',
         r'\label{tab:mpsnr}', r'\centering', r'\scriptsize', r'\setlength{\tabcolsep}{2.6pt}',
         r'\begin{tabular}{c' + 'c' * len(M7) + '}', r'\toprule', 'Case & ' + ' & '.join(SHORT[m] for m in M7) + r'\\']
for c in CUBES:
    lines += [r'\midrule', r'\multicolumn{%d}{l}{\emph{%s} (%s)}\\' % (len(M7) + 1, CN[c], SIZE[c])]
    for case in range(1, 11):
        full = cnt.get((c, case, 'SC'), 0)
        v = [cm.get((m, c, case), np.nan) if (m == 'FastHyMix' or cnt.get((c, case, m), 0) == full) else np.nan for m in M7]
        b = bold_best(v)
        b = [(r'$\times$' if m == 'FastHyMix' and cnt.get((c, case, m), 0) == 0 and 'FastHyMix' in set(df.meth) else x)
             + (DAG if m == 'FastHyMix' and 0 < cnt.get((c, case, m), 0) < full else '') for x, m in zip(b, M7)]
        lines.append(' & '.join([str(case)] + b) + r'\\')
lines += [r'\bottomrule', r'\end{tabular}'] + (note(DAGNOTE_CASE) if FHFAIL else []) + [r'\end{table}']
write('tab_mpsnr.tex', lines)

# ------------------------------------------------------------------ averages over all cells + time (one column)
s0df = df[(df.seed == 0) & (df.case <= 8)]
lines = [r'\begin{table}[!t]',
         r'\caption{Mean over all cases and realizations of the simulated cubes. Time: one run on case~8 (seed~0) on an otherwise idle machine. Below the inner rule: S3TTV, run on seed~0, cases 1--8 only (about ' + S3MIN + r'\,min per cube), and SC-RCCTV on the same eight inputs.}',
         r'\label{tab:avg}', r'\centering', r'\scriptsize', r'\setlength{\tabcolsep}{2.0pt}', r'\begin{tabular}{llccccc}', r'\toprule',
         r'Cube & Method & MPSNR$\uparrow$ & MSSIM$\uparrow$ & SAM$\downarrow$ & ERGAS$\downarrow$ & Time (s)\\']
for c in CUBES:
    sub = df[(df.cube == c) & df.meth.isin(M7)]
    g = sub.groupby('meth')[['mpsnr', 'mssim', 'sam', 'ergas']].mean()
    done = [m for m in M7 if m in g.index and complete(m, c)]
    gb = g.loc[done]
    lines.append(r'\midrule')
    for i, m in enumerate(M7):
        t = tsec(c + '_c8', m)
        ts = TBD if np.isnan(t) else fsec(t)
        lead = CN[c].split()[0] if i == 0 else ''
        if m not in done:
            lines.append(f"{lead} & {NAME[m]} & {' & '.join([TBD] * 4)} & {ts}" + r'\\')
            continue
        r = g.loc[m]
        f = lambda x, b, fm: (r'\textbf{' + fm.format(x) + '}') if x == b else fm.format(x)
        dag = DAG if len(sub[sub.meth == m]) < int(ncell[c]) else ''
        lines.append(f"{lead} & {NAME[m]}{dag} & {f(r.mpsnr, gb.mpsnr.max(), '{:.2f}')} & {f(r.mssim, gb.mssim.max(), '{:.4f}')} & "
                     f"{f(r.sam, gb.sam.min(), '{:.2f}')} & {f(r.ergas, gb.ergas.min(), '{:.2f}')} & {ts}" + r'\\')
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
        lines.append(f" & {NAME[m]}$^\\ddagger$ & {vals} & {ts}" + r'\\')
lines += [r'\bottomrule', r'\end{tabular}'] + note((DAGNOTE + ' ' if FHFAIL else '') + r'$^\ddagger$Seed~0, cases 1--8.') + [r'\end{table}']
write('tab_avg.tex', lines)

# ------------------------------------------------------------------ GRSL: compact table by noise family (two columns)
lines = [r'\begin{table*}[!t]',
         r'\caption{Simulated data: MPSNR (dB) by noise family (G+I: Gaussian and impulse noise, cases 1--3; Str.: column stripes, cases 4--8; B.-v.: band-varying noise, cases 9--10) and over all cells, MSSIM and SAM over all cells; mean over six realizations. Time: case~8, Pavia University/Beltsville. Best per column in bold. Per-case values: Table~S1 of the supplementary material.}',
         r'\label{tab:family}', r'\centering', r'\scriptsize', r'\setlength{\tabcolsep}{3.1pt}',
         r'\begin{tabular}{l' + 'cccccc' * 2 + 'c}', r'\toprule',
         r' & \multicolumn{6}{c}{Pavia University (' + SIZE['pavia120'] + r')} & \multicolumn{6}{c}{Beltsville (' + SIZE['belt100'] + r')} & \\',
         r'\cmidrule(lr){2-7}\cmidrule(lr){8-13}',
         r'Method & ' + ' & '.join(['G+I', 'Str.', 'B.-v.', 'All', 'MSSIM', 'SAM'] * 2) + r' & Time (s)\\', r'\midrule']
vals = {}
for c in CUBES:
    sub = df[(df.cube == c) & df.meth.isin(M7)]
    for m in M7:
        if not complete(m, c):
            vals[(c, m)] = [np.nan] * 6
            continue
        s = sub[sub.meth == m]
        vals[(c, m)] = [s[s.case.isin(cs)].mpsnr.mean() for cs in FAM.values()] + [s.mpsnr.mean(), s.mssim.mean(), s.sam.mean()]
best = {}
for c in CUBES:
    for j in range(6):
        col = np.array([vals[(c, m)][j] for m in M7], float)
        best[(c, j)] = np.nanmin(col) if j == 5 else np.nanmax(col)
for m in M7:
    cells = []
    for c in CUBES:
        for j, x in enumerate(vals[(c, m)]):
            if np.isnan(x):
                cells.append(TBD)
                continue
            s = f'{x:.4f}' if j == 4 else f'{x:.2f}'
            s = r'\textbf{' + s + '}' if x == best[(c, j)] else s
            fam_hit = j >= 3 or bool(FAILCASES[c] & set(list(FAM.values())[j]))   # family column touched by a crash; All/MSSIM/SAM always
            if m == 'FastHyMix' and FAILCASES[c] and fam_hit:
                s += DAG
            cells.append(s)
    t = [tsec(c + '_c8', m) for c in CUBES]
    cells.append(TBD if any(np.isnan(t)) else '/'.join(fsec(x) for x in t))
    lines.append(NAME[m] + ' & ' + ' & '.join(cells) + r'\\')
lines += [r'\bottomrule', r'\end{tabular}'] + (note(DAGNOTE) if FHFAIL else []) + [r'\end{table*}']
write('tab_family.tex', lines, [GRSL])

# ------------------------------------------------------------------ seed 0, cases 1-8 with S3TTV (two columns)
m3 = ['E3DTV', 'FastHyMix', 'tCTVd', 'S3TTV', 'SC']
lines = [r'\begin{table*}[!t]',
         r'\caption{MPSNR (dB) on the official inputs (seed~0, cases 1--8), the only cells on which S3TTV was run. $\times$: the released code stopped with an error.}',
         r'\label{tab:seed0}', r'\centering', r'\footnotesize', r'\setlength{\tabcolsep}{4pt}',
         r'\begin{tabular}{c' + 'c' * len(m3) + 'c' * len(m3) + '}', r'\toprule',
         r' & \multicolumn{%d}{c}{Pavia University} & \multicolumn{%d}{c}{Beltsville}\\' % (len(m3), len(m3)),
         r'\cmidrule(lr){2-%d}\cmidrule(lr){%d-%d}' % (len(m3) + 1, len(m3) + 2, 2 * len(m3) + 1),
         'Case & ' + ' & '.join(NAME[m] for m in m3) + ' & ' + ' & '.join(NAME[m] for m in m3) + r'\\', r'\midrule']
ran = set(df.meth)
s0 = _ns['s0']
for case in range(1, 9):
    row = [str(case)]
    for c in CUBES:
        v = [pv.loc[(c, case, 0), m] if (m in pv and (c, case, 0) in pv.index) else np.nan for m in m3]
        b = bold_best(v)
        row += [(r'$\times$' if (np.isnan(val) and m in ran) else x) for x, m, val in zip(b, m3, v)]
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
lines += [r'\bottomrule', r'\end{tabular}', r'\end{table*}']
write('tab_seed0.tex', lines)

# ------------------------------------------------------------------ ablation (one column)
fam = {'G+I': [1, 2, 3], 'Str.': [4, 5, 6, 7, 8], 'B.-v.': [9, 10]}
variants = [('SC', 'SC-RCCTV (full)'), ('SC_noT', 'without stripe term $T$'), ('SC_noG', 'without dense term $G$'),
            ('SC_scalar', r'scalar $\sigma$ instead of $\sigma_k$'), ('SC_fix', r'$\hat\sigma_k^{(0)}$ kept fixed (two-stage)'),
            ('SCv', 'column stripes only'), ('SC_true', r'true $\sigma_k$ (all else unchanged)'),
            ('OLD', r'RCCTV-K, true $\sigma$ and radius'), ('OLDband', r'RCCTV-K, true $\sigma_k$ and radius')]
cap = (r'Ablation: mean MPSNR gain (dB) over RCCTV by noise family (G+I: cases 1--3; Str.: 4--8; B.-v.: 9--10; both cubes, all realizations). '
       r'Below the rule: SC-RCCTV with the true $\sigma_k$ in place of the self-calibrated ones, all else unchanged; and RCCTV-K, a different noise-aware model: RCCTV with a dense term weighted by the true noise level (scalar $\sigma$, i.e.\ the root mean square of $\sigma_k$ in cases 9--10, or the band-wise $\sigma_k$) and an $\ell_1$-ball stripe set with the true radius.')
for fname, keep, dirs in (('tab_abl.tex', ['SC', 'SC_noT', 'SC_noG', 'SC_scalar', 'SC_fix', 'SC_true', 'OLD', 'OLDband'], [GRSL]),
                          ('tab_abl_full.tex', [k for k, _ in variants], BOTH)):
    lines = [r'\begin{table}[!t]', r'\caption{' + cap + '}', r'\label{tab:abl}' if fname == 'tab_abl.tex' else r'\label{tab:ablfull}',
             r'\centering', r'\scriptsize', r'\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{lcccc}', r'\toprule',
             r'Variant & ' + ' & '.join(fam) + r' & All\\', r'\midrule']
    for key, lab in variants:
        if key not in keep:
            continue
        if key == 'SC_true':
            lines.append(r'\midrule')
        if key not in pv:
            lines.append(lab + ' & ' + ' & '.join([TBD] * 4) + r'\\')
            continue
        d = (pv[key] - pv['RCCTV']).dropna()
        cells = [sgn(d[d.index.get_level_values('case').isin(cs)].mean()) for cs in fam.values()] + [sgn(d.mean())]
        lines.append(lab + ' & ' + ' & '.join(cells) + r'\\')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    write(fname, lines, dirs)

# ------------------------------------------------------------------ real data (one column)
order = ['Input', 'RCCTV', 'LRTDTV', 'E3DTV', 'tCTVd', 'FastHyMix', 'L1HyMixDe', 'SC']
lines = [r'\begin{table}[!t]',
         r'\caption{Real data. Indian Pines: overall accuracy (OA, mean $\pm$ std over 5 splits) and $\kappa$ of an RBF-SVM with 10\% training samples per class. Urban ($200{\times}200$): stripe index averaged over the 20 noisiest bands of the input, and largest change of a band mean, both $\times10^{-3}$. Time: Indian Pines, idle machine; t-CTV was run on Indian Pines only.}',
         r'\label{tab:real}', r'\centering', r'\scriptsize', r'\setlength{\tabcolsep}{3pt}', r'\begin{tabular}{lccccc}', r'\toprule',
         r'Method & OA (\%)$\uparrow$ & $\kappa\uparrow$ & Stripe$\downarrow$ & Shift$\downarrow$ & Time (s)\\', r'\midrule']
for k in order:
    u = rm[(rm.ds == 'urban200') & (rm.meth == k)]
    oa = f"{cl.loc[k, 'OA']:.1f}\\,$\\pm$\\,{cl.loc[k, 'OA_std']:.1f}" if k in cl.index else '--'
    kp = f"{cl.loc[k, 'kappa']:.3f}" if k in cl.index else '--'
    si = f"{1e3 * u.stripeIdx_top20.iloc[0]:.1f}" if len(u) else '--'
    sh = ('0' if k == 'Input' else f"{1e3 * u.shift_max.iloc[0]:.1f}") if len(u) else '--'
    t = tsec('indian220', k)
    sec = '--' if k == 'Input' else (TBD if np.isnan(t) else fsec(t))
    lines.append(f'{NAME[k]} & {oa} & {kp} & {si} & {sh} & {sec}' + r'\\')
lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
write('tab_real.tex', lines)

# ------------------------------------------------------------------ full-size scenes (one column)
LM = ['RCCTV', 'LRTDTV', 'E3DTV', 'tCTVd', 'FastHyMix', 'L1HyMixDe', 'SC']
lines = [r'\begin{table}[!t]',
         r'\caption{Full-size scenes, one run per method on an otherwise idle machine: Pavia University $610{\times}340{\times}103$ (cases 8 and 10) and Urban $307{\times}307{\times}210$ (real; stripe index averaged over the 20 bands with the highest stripe index in the input, $\times10^{-3}$). S3TTV is not included: its memory use exceeds 16\,GB at this size.}',
         r'\label{tab:large}', r'\centering', r'\scriptsize', r'\setlength{\tabcolsep}{3pt}', r'\begin{tabular}{lccccc}', r'\toprule',
         r' & \multicolumn{3}{c}{Pavia University} & \multicolumn{2}{c}{Urban}\\', r'\cmidrule(lr){2-4}\cmidrule(lr){5-6}',
         r'Method & c8 & c10 & Time (s) & Stripe$\downarrow$ & Time (s)\\', r'\midrule']
p8 = [lget('paviaU_c8', m, 'mpsnr') for m in LM]
p10 = [lget('paviaU_c10', m, 'mpsnr') for m in LM]
ust = [lget('urban307', m, 'stripe_top20') * 1e3 for m in LM]
b8, b10, bst = bold_best(p8), bold_best(p10), bold_best(ust, higher=False, fmt='{:.1f}')
for i, m in enumerate(LM):
    tp = np.nanmean([lget('paviaU_c8', m, 'sec'), lget('paviaU_c10', m, 'sec')])
    tu = lget('urban307', m, 'sec')
    lines.append(f"{NAME[m]} & {b8[i]} & {b10[i]} & {TBD if np.isnan(tp) else fsec(tp)} & {bst[i]} & {TBD if np.isnan(tu) else fsec(tu)}" + r'\\')
lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
write('tab_large.tex', lines)

for d in BOTH:
    shutil.copy('paper_rsl/numbers.tex', d + 'numbers.tex')
print('IEEE tables written to', ', '.join(BOTH))
