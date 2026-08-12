#!/usr/bin/env python3
"""Generate the two paper figures in base-LaTeX picture-environment code."""
import math

HOLDS = {2: [1, 2, 4], 3: [1, 3], 4: [1, 2, 4], 5: [1],
         6: [1, 2, 3, 4, 6], 7: [1], 8: [1, 2, 4], 9: [1, 3],
         10: [1, 2, 4], 11: [1], 12: [1, 2, 3, 4, 6]}
MMIN, MMAX, FMAX = 2, 12, 8
ROWS = MMAX - MMIN + 1          # 11
Y = lambda m: MMAX - m          # row m -> y of its cell bottom


def fig_grid():
    L = []
    A = L.append
    A(r'\setlength{\unitlength}{0.66cm}%')
    A(r'\begin{picture}(10.1,13.3)(-2.05,-0.1)')
    A(r'  % --- shaded cells: sieving holds')
    for m in range(MMIN, MMAX + 1):
        for f in HOLDS[m]:
            A(r'  \put(%d,%d){\textcolor{cellpale}{\rule{1\unitlength}{1\unitlength}}}'
              % (f - 1, Y(m)))
    A(r'  % --- grid')
    A(r'  \thinlines')
    for j in range(ROWS + 1):
        A(r'  \put(0,%d){\textcolor{gridgray}{\line(1,0){%d}}}' % (j, FMAX))
    for i in range(FMAX + 1):
        A(r'  \put(%d,0){\textcolor{gridgray}{\line(0,1){%d}}}' % (i, ROWS))
    A(r'  % --- divider: known cases (f=1,2) | new cases')
    A(r'  \put(2,0){\textcolor{rulegray}{\rule{0.9pt}{%d\unitlength}}}' % ROWS)
    A(r'  % --- verdict marks')
    for m in range(MMIN, MMAX + 1):
        for f in HOLDS[m]:
            A(r'  \put(%d,%d){\makebox(1,1){\textcolor{azurite}{$\bullet$}}}'
              % (f - 1, Y(m)))
    A(r'  % --- axis ticks')
    for i in range(FMAX):
        A(r'  \put(%d,%s){\makebox(1,0.55)[b]{\footnotesize $%d$}}'
          % (i, '11.13', i + 1))
    for m in range(MMIN, MMAX + 1):
        A(r'  \put(-1.42,%d){\makebox(1.25,1)[r]{\footnotesize $%d$}}' % (Y(m), m))
    A(r'  % --- axis titles')
    A(r'  \put(0,12.34){\makebox(%d,0.6)[b]{\small fixed-point count $f$}}' % FMAX)
    A(r'  \put(-2.05,0){\makebox(0.6,%d){\rotatebox{90}{\small order $m$}}}' % ROWS)
    A(r'\end{picture}')
    return '\n'.join(L)


ORDERS = [(1, 'allowed'), (2, 'allowed'), (3, 'allowed'), (4, 'allowed'),
          (5, 'forbidden'), (6, 'allowed'), (7, 'forbidden'), (8, 'forbidden')]
STEP, RAD = 1.70, 0.62


def num(x):
    """Fixed-point decimal: picture coordinates must never be in E notation."""
    s = '%.4f' % x
    s = s.rstrip('0').rstrip('.')
    return '0' if s in ('', '-0', '-') else s


def fig_orders():
    L = []
    A = L.append
    A(r'\setlength{\unitlength}{1cm}%')
    A(r'\begin{picture}(13.2,2.15)(-0.66,-1.48)')
    for idx, (f, status) in enumerate(ORDERS):
        cx = idx * STEP
        col = 'azurite' if status == 'allowed' else 'oxide'
        A(r'  %% --- f = %d (%s)' % (f, status))
        A(r'  \put(%s,0){\textcolor{gridgray}{\circle{%s}}}' % (num(cx), num(2 * RAD)))
        for i in range(f):
            th = math.radians(90 + 360.0 * i / f)
            A(r'  \put(%s,%s){\textcolor{%s}{\circle*{0.16}}}'
              % (num(cx + RAD * math.cos(th)), num(RAD * math.sin(th)), col))
        A(r'  \put(%s,-1.04){\makebox(1,0.42)[b]{\footnotesize $%d$}}' % (num(cx - 0.5), f))
        A(r'  \put(%s,-1.44){\makebox(1.6,0.36)[b]{\textcolor{%s}{\tiny %s}}}'
          % (num(cx - 0.8), col, status))
    A(r'\end{picture}')
    return '\n'.join(L)


if __name__ == '__main__':
    print('%%%% FIG-GRID %%%%')
    print(fig_grid())
    print('%%%% FIG-ORDERS %%%%')
    print(fig_orders())
