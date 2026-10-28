"""SIDE reading for ibm_marrakesh (pairs_repeat/PREREG_PAIRS.md): the same R1-R5, by control_run/analyze_pairs_iq.py UNCHANGED, on the pairs
NOT wired (cz, the run's own saved properties) to measured qubits 82, 100, 119, whose 0/1 IQ clusters overlap or miss on marrakesh
(flipcheck_run/derive_flip_guard.txt: flip-rule failure 98.9 %, 76.8 %, 71.4 % on 21 Oct calibration-1 data). Printed BESIDE the main
reading; it does not replace it. Heuremen LLC, Wayfinder and Bones, 22 October 2026.
How: copies the run's marrakesh job records into a temp folder with the excluded pairs removed from 'pairs' and their two fin columns removed
from every fin_ register (kept pairs keep their order, singles follow, as the runner lays them out); every other register (c_, m<k>_) and the
properties file are copied as they are. Then runs the unchanged analysis on that folder. R5 still counts all 39 measured qubits (82, 100,
119 included) because the analysis computes it over E; read R5 from the main reading.
Usage: side_marrakesh.py [FOLDER] [TAG]   (default FOLDER = this folder; TAG names the output side_marrakesh_<TAG>.txt)"""
import os, sys, json, subprocess, tempfile, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ANALYSIS = os.path.join(os.path.dirname(HERE), 'control_run', 'analyze_pairs_iq.py')
SRC = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE; TAG = sys.argv[2] if len(sys.argv) > 2 else 'repeat'
CHIP = 'ibm_marrakesh'; BAD = {82, 100, 119}
props_f = os.path.join(SRC, 'pairs_backend_properties_%s.json' % CHIP); props = json.load(open(props_f))
adj = {}
for g in props['gates']:
    if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
tmp = tempfile.mkdtemp(prefix='side_marrakesh_')
try:
    shutil.copy(props_f, tmp); lines = []
    for c in (0, 1):
        jf = os.path.join(SRC, 'pairsjob_%s_copy%d.json' % (CHIP, c)); rf = os.path.join(SRC, 'pairsraw_%s_copy%d.npz' % (CHIP, c))
        if not (os.path.exists(jf) and os.path.exists(rf)): lines.append('copy %d: not on disk' % c); continue
        rec = json.load(open(jf)); P = rec['pairs']; nP = len(P); nS = len(rec['singles'])
        drop = [j for j, (a, b) in enumerate(P) if (adj.get(a, set()) | adj.get(b, set())) & BAD]
        keepj = [j for j in range(nP) if j not in set(drop)]
        cols = [x for j in keepj for x in (2 * j, 2 * j + 1)] + list(range(2 * nP, 2 * nP + nS))
        lines.append('copy %d: %d pairs, %d removed (wired to %s): %s' % (c, nP, len(drop), sorted(BAD), [P[j] for j in drop]))
        z = np.load(rf); out = {}
        for k in z.files:
            a = z[k]
            if k.startswith('fin_'): assert a.shape[1] == 2 * nP + nS; a = a[:, cols]
            out[k] = a
        json.dump(dict(rec, pairs=[P[j] for j in keepj]), open(os.path.join(tmp, os.path.basename(jf)), 'w'))
        np.savez(os.path.join(tmp, os.path.basename(rf)), **out)
    res = subprocess.run([sys.executable, ANALYSIS, tmp], capture_output=True, text=True)
    txt = ('SIDE reading, ibm_marrakesh, pairs not wired to measured qubits 82/100/119 (side_marrakesh.py; analysis = control_run/analyze_pairs_iq.py '
           'unchanged).\nSource: %s\n%s\nR5 below still includes 82/100/119 (computed over all 39 measured qubits).\n' % (SRC, '\n'.join(lines)) + res.stdout + res.stderr)
    open(os.path.join(HERE, 'side_marrakesh_%s.txt' % TAG), 'w').write(txt); print(txt)
finally: shutil.rmtree(tmp)
