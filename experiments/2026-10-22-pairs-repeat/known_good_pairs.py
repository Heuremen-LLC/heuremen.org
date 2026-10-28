"""KNOWN-GOOD check for the pairs repeat (pairs_repeat/PREREG_PAIRS.md): the idle arm at N = 0 (no readouts, no idle steps) must show the Bell
pairs entangled. analyze_pairs_iq.py does not print the N = 0 arm, so this small reader does, with the analysis's own classification
(mu0, mu1 = mean of all calibration-0 / -1 shots, both copies; bit = 1 if the projection on mu1 - mu0 exceeds 0.5) and its concurrence
formula copied verbatim. Wired = at least one measured qubit cz-coupled to either member (the run's saved properties), as in R1.
Read with the analysis's own ruler: median N = 0 concurrence over wired pairs, 95 % bootstrap (2000), four outcomes; known-good holds if ABOVE 0.
Heuremen LLC, Wayfinder and Bones, 22 October 2026.   Usage: known_good_pairs.py [FOLDER] [TAG] -> known_good_<TAG>.txt here"""
import os, sys, json, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE; TAG = sys.argv[2] if len(sys.argv) > 2 else 'repeat'
SETS = ('ZZ', 'XX', 'YY', 'XY', 'YX')
def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))
out = ['KNOWN-GOOD (N = 0 idle concurrence on wired pairs), source %s' % SRC]
chips = sorted({os.path.basename(f).split('_copy')[0][len('pairsjob_'):] for f in glob.glob(os.path.join(SRC, 'pairsjob_*_copy*.json'))})
for chip in chips:
    copies = {}
    for c in (0, 1):
        jf = os.path.join(SRC, 'pairsjob_%s_copy%d.json' % (chip, c)); rf = os.path.join(SRC, 'pairsraw_%s_copy%d.npz' % (chip, c))
        if os.path.exists(jf) and os.path.exists(rf): copies[c] = (json.load(open(jf)), np.load(rf))
    if not copies: out.append('%s NOTHING_RETURNED' % chip); continue
    rec0 = next(iter(copies.values()))[0]; E = set(rec0['E']); P = rec0['pairs']
    props = json.load(open(os.path.join(SRC, 'pairs_backend_properties_%s.json' % chip))); adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    wired = np.array([bool((adj.get(a, set()) | adj.get(b, set())) & E) for a, b in P])
    cal = {s: np.concatenate([z['c_%d' % i] for (r, z) in copies.values() for i, m in enumerate(r['meta']) if m['arm'] == 'CAL' and m['state'] == s]) for s in '01'}
    mu0 = cal['0'].mean(0); u = cal['1'].mean(0) - mu0; proj = lambda x, q: np.real((x - mu0[q]) * np.conj(u[q])) / np.abs(u[q]) ** 2
    for c, (r, z) in sorted(copies.items()):
        cells = [dict() for _ in P]
        for i, m in enumerate(r['meta']):
            if not (m['arm'] == 'P' and m['mode'] == 'D' and m['N'] == 0): continue
            fin = z['fin_%d' % i]; st = m['setting']
            for j, (a, b) in enumerate(P):
                b1 = (proj(fin[:, 2 * j], a) > 0.5).astype(int); b2 = (proj(fin[:, 2 * j + 1], b) > 0.5).astype(int)
                cells[j][st] = [np.mean((b1 == x) & (b2 == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))] if st == 'ZZ' else np.mean((1 - 2 * b1) * (1 - 2 * b2))
        cv = np.array([conc(x) if len(x) == 5 else np.nan for x in cells])[wired]
        x = cv[np.isfinite(cv)]; rng = np.random.default_rng(20261028)   # the analysis's own ruler: median, 95 % bootstrap (2000), four outcomes
        med = float(np.median(x)); b = np.median(rng.choice(x, (2000, len(x))), 1); lo, hi = np.percentile(b, [2.5, 97.5])
        st = 'CANNOT_TELL' if len(x) < 3 else 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
        out.append('  %-14s copy %d: wired pairs %d, N=0 idle concurrence median %.3f [%.3f, %.3f] -> %s'
                   % (chip, c, len(x), med, lo, hi, st))
txt = '\n'.join(out); open(os.path.join(HERE, 'known_good_%s.txt' % TAG), 'w').write(txt + '\n'); print(txt)
