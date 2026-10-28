"""Synthetic test of analyze_flip.py in the real raw-IQ schema (PLAN_FLIP.md). Heuremen LLC, Wayfinder and Bones, 22 October 2026.
Roles (E, S), cz wiring and every qubit's 0/1 IQ cluster centres and spreads are the REAL ibm_kingston ones from the 21 Oct raw-readout run
(control_run/iqjob_ibm_kingston.json, iqraw_ibm_kingston.npz, iq_backend_properties_ibm_kingston.json). Same pubs, meta, register names,
shapes and dtype as run_flip.py writes (26 pubs, 3000 shots). Planted:
  (a) loss 0.008 on every untouched qubit in the readout arm, no hidden leak
  (b) loss caused ONLY by hidden-leaked measured neighbours: each measured qubit leaks at its first mid readout in 1 % of readout-arm shots,
      its IQ point then sits inside its own 0 cluster on every later read (and after the X: X does not move it), and every untouched qubit
      wired to a leaked neighbour loses 0.8 coherence in that shot
  (c) null: nothing planted; six seeds
Common to all: idle coherence 0.90-0.98 per qubit times exp(-0.01 N), a fixed phase per qubit, 0.5 % visible readout excitation of measured
qubits (lands in the 1 cluster; the selection sees it), 1.5 % baseline flip failure in BOTH arms (T1 during the final read; lands in 0).
Usage: synth_flip.py [scenario ...]   (default: a b b1..b5 null0..null5; bK = scenario b, seed offset K). Writes into a temp folder, runs analyze_flip.py, deletes the raw file."""
import os, sys, json, subprocess, tempfile, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CTRL = os.path.join(os.path.dirname(HERE), 'control_run'); CHIP = 'ibm_kingston'
SHOTS = 3000; SEED = 20261023; NS = (4, 8)
if '--chip' in sys.argv: i = sys.argv.index('--chip'); CHIP = sys.argv[i + 1]; del sys.argv[i:i + 2]
rec0 = json.load(open(os.path.join(CTRL, 'iqjob_%s.json' % CHIP))); E = rec0['E']; S = rec0['S']; n = len(E) + len(S)
props = json.load(open(os.path.join(CTRL, 'iq_backend_properties_%s.json' % CHIP))); cz = [g for g in props['gates'] if g['gate'] == 'cz']
zr = np.load(os.path.join(CTRL, 'iqraw_%s.npz' % CHIP))
C0 = np.concatenate([zr['c_%d' % i] for i, m in enumerate(rec0['meta']) if m['arm'] == 'CAL' and m['state'] == '0'])
C1 = np.concatenate([zr['c_%d' % i] for i, m in enumerate(rec0['meta']) if m['arm'] == 'CAL' and m['state'] == '1'])
MU0 = np.median(C0.real, 0) + 1j * np.median(C0.imag, 0); MU1 = np.median(C1.real, 0) + 1j * np.median(C1.imag, 0)
SIG = np.median(np.abs(C0 - MU0), 0) / np.sqrt(2 * np.log(2))      # per-quadrature spread of the 0 cluster (median radius of a 2-D Gaussian = 1.1774 sigma)
adj = {}
for g in cz: a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
Epos = {e: i for i, e in enumerate(E)}; NB = [[Epos[e] for e in adj.get(s, ()) if e in Epos] for s in S]

def pts(rng, qubits, lev):   # lev: shots x len(qubits) in {0,1,2}; 2 = leaked, placed 1 sigma off the 0 centre (inside the 0 cluster)
    q = np.asarray(qubits); mu = np.where(lev == 0, MU0[q], np.where(lev == 1, MU1[q], MU0[q] + SIG[q] * (1 + 0j)))
    return (mu + SIG[q] * (rng.standard_normal(lev.shape) + 1j * rng.standard_normal(lev.shape))).astype(np.complex64)

def make(scn, seed, out):
    rng = np.random.default_rng(seed); g = np.random.default_rng(SEED)        # g: fixed qubit properties
    c0 = g.uniform(0.90, 0.98, len(S)); ph = g.uniform(0, 2 * np.pi, len(S))
    meta, data = [], {}
    base = [dict(arm='CAL', state='0'), dict(arm='CAL', state='1'), dict(arm='CAL', state='0X'), dict(arm='S', mode='D', N=0, setting='XX'), dict(arm='S', mode='D', N=0, setting='YY')]
    base += [dict(arm='S', mode=mo, N=N, setting=st) for N in NS for mo in ('D', 'M') for st in ('XX', 'YY')]
    meta = [dict(m, copy=c) for c in (0, 1) for m in base]; order = list(range(len(meta))); rng.shuffle(order); meta = [meta[i] for i in order]
    Q = np.arange(n); nE = len(E)
    for i, m in enumerate(meta):
        if m['arm'] == 'CAL':
            if m['state'] in ('0', '1'):
                lev = np.full((SHOTS, n), int(m['state']))
                if m['state'] == '0': lev[rng.random((SHOTS, n)) < 0.005] = 1
                else: lev[rng.random((SHOTS, n)) < 0.015] = 0
                data['c_%d' % i] = pts(rng, Q, lev)
            else:
                l0 = np.zeros((SHOTS, n), int); l0[rng.random((SHOTS, n)) < 0.005] = 1; l1 = 1 - l0; l1[rng.random((SHOTS, n)) < 0.015] = 0
                data['c0_%d' % i] = pts(rng, Q, l0); data['c1_%d' % i] = pts(rng, Q, l1)
            continue
        N, mo = m['N'], m['mode']
        el = np.zeros((SHOTS, nE), int)                                           # measured qubits' level during the steps
        if mo == 'M':
            el[rng.random((SHOTS, nE)) < 0.005] = 1                              # visible readout excitation
            if scn == 'b': el[rng.random((SHOTS, nE)) < 0.01] = 2                 # hidden leak at the first readout
        for k in range(N):
            data['m%d_%d' % (k, i)] = pts(rng, E, el) if mo == 'M' else np.zeros((SHOTS, nE), np.complex64)
        ef = np.where(el == 0, 1, np.where(el == 1, 0, 2)); ef[(ef == 1) & (rng.random((SHOTS, nE)) < 0.015)] = 0   # X, then T1 in the final read
        data['efin_%d' % i] = pts(rng, E, ef)
        c = np.broadcast_to(c0 * np.exp(-0.01 * N), (SHOTS, len(S))).copy()
        if mo == 'M' and scn == 'a': c -= 0.008
        if mo == 'M' and scn == 'b':
            leaked = el == 2; hit = np.stack([leaked[:, x].any(1) if x else np.zeros(SHOTS, bool) for x in NB], 1); c[hit] -= 0.8
        cosv = np.cos(ph) if m['setting'] == 'XX' else np.sin(ph)
        bit = (rng.random((SHOTS, len(S))) > (1 + c * cosv) / 2).astype(int)
        data['fin_%d' % i] = pts(rng, S, bit)
    json.dump({'chip': CHIP, 'job_id': 'SYNTH-%s-%d' % (scn, seed), 'meta': meta, 'order': order, 'E': E, 'S': S, 'shots': SHOTS, 'meas_type': 'kerneled'},
              open(os.path.join(out, 'flipjob_%s.json' % CHIP), 'w'))
    json.dump({'gates': cz}, open(os.path.join(out, 'flip_backend_properties_%s.json' % CHIP), 'w'))
    np.savez(os.path.join(out, 'flipraw_%s.npz' % CHIP), **data)

if __name__ == '__main__':
    scns = sys.argv[1:] or ['a', 'b'] + ['b%d' % k for k in range(1, 6)] + ['null%d' % k for k in range(6)]
    for sc in scns:
        out = tempfile.mkdtemp(prefix='synthflip_'); seed = SEED + 100 + (int(sc[4:]) if sc.startswith('null') else ord(sc[0]) + 1000 * int(sc[1:] or 0))
        make('null' if sc.startswith('null') else sc[0], seed, out)
        r = subprocess.run([sys.executable, os.path.join(HERE, 'analyze_flip.py'), out], capture_output=True, text=True)
        print('#' * 30, 'SCENARIO', sc, 'seed', seed); print(r.stdout + r.stderr, flush=True); shutil.rmtree(out)
