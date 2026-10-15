"""Create experiments/matlab/tctv_trpca_lam.m from the official t-CTV code (not redistributed here).

The variant differs from the official TCTV_TRPCA.m in exactly two places: the function name, and one added line
that scales the default sparse weight by opts.lambda_mult (selected on development cubes; see tune_tctv.m).
Usage:  python tools/make_tctv_trpca_lam.py code/official_tctv/TCTV_TRPCA.m
"""
import os
import sys

src = sys.argv[1] if len(sys.argv) > 1 else os.path.join('code', 'official_tctv', 'TCTV_TRPCA.m')
dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments', 'matlab', 'tctv_trpca_lam.m')
text = open(src, encoding='utf-8', errors='replace').read()
head = 'function [X, E, obj, err, iter] = TCTV_TRPCA(M, opts)'
anchor = "if isfield(opts, 'detail');"
if head not in text or anchor not in text:
    sys.exit('unexpected TCTV_TRPCA.m: function header or opts.detail line not found')
text = text.replace(head, 'function [X, E, obj, err, iter] = tctv_trpca_lam(M, opts)', 1)
lines = text.splitlines(keepends=True)
i = next(k for k, l in enumerate(lines) if l.lstrip().startswith(anchor))
nl = '\r\n' if lines[i].endswith('\r\n') else '\n'
lines.insert(i + 1, "if isfield(opts, 'lambda_mult');      lambda             = lambda*opts.lambda_mult; end  % ADDED (only change): scale default lambda" + nl)
lines.insert(0, '% Copy of official TCTV_TRPCA.m (code/official_tctv) with ONE added line: opts.lambda_mult scales the default lambda.' + nl)
open(dst, 'w', encoding='utf-8', newline='').write(''.join(lines))
print('written', os.path.normpath(dst))
