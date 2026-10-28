"""Edges IBM states as zz = 0.00 in the flat-beat job (ibm_phoenix): the committed fit's own grid and envelope applied to the three
edges the committed output did not fit. Defined after the data (an outside referee asked for them as an extra line), per copy and pooled.
Run from this folder: python3 zero_edges_postdata.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '2026-10-22-one-neighbour-revival'))   # analyze_zz, as analyze_flat's '../zz_revival'
sys.argv = [sys.argv[0], HERE]
import numpy as np, analyze_flat as A
Z = A.Z; J = A.J2('ibm_phoenix')
T = np.array(J.taus) * 1e-6; floor = 5 * np.sqrt(2.0 / (np.pi * 2 * J.rec['shots']))
zz = Z.zz_table(J.props, J.edges)
def fit(arm, j, copies):
    Cz = np.array([Z.coh(J, 'Z', t, copies)[0] for t in J.taus])[:, j]; Cp = np.array([Z.coh(J, arm, t, copies)[0] for t in J.taus])[:, j]
    lv = Cz > floor; m = Cz[lv][None, :] * np.abs(np.cos(np.pi * A.F[:, None] * T[lv][None, :]))
    return A.F[np.argmin(((m - Cp[lv][None, :]) ** 2).sum(1))] / 1e3, int(lv.sum())
for j, (b1, a, b2) in enumerate(J.trip):
    for b, arm in ((b1, 'P1'), (b2, 'P2')):
        if zz.get((min(a, b), max(a, b)), np.nan) == 0:
            f, n = fit(arm, j, (0, 1)); f0, _ = fit(arm, j, (0,)); f1, _ = fit(arm, j, (1,))
            print(f'edge {a}-{b}: stated zz 0.00 kHz; fitted beat {f:.2f} kHz (copy 0 {f0:.2f}, copy 1 {f1:.2f}; {n} live delays)')
