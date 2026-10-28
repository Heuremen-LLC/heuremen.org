"""Synthetic test of analyze.py: builds every run's real circuits list on FakeFez (offline), plants known effects, writes the same files a
fired job writes, runs the analysis. Planted values are printed first so the read-back can be compared."""
import os, sys, json, random, warnings; warnings.filterwarnings('ignore'); sys.argv = sys.argv[:1]
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import runs, analyze
from qiskit_ibm_runtime.fake_provider import FakeFez
NULL = os.environ.get('NULL') == '1'; SEEDN = int(os.environ.get('SEED', '5'))
b = FakeFez(); n = b.num_qubits; rng = np.random.default_rng(SEEDN); OUT = os.path.join(HERE, 'synth_null_%d' % SEEDN if NULL else 'synth'); os.makedirs(OUT, exist_ok=True)
analyze.RNG = np.random.default_rng(SEEDN + 1000)
props = b.properties().to_dict(); json.dump(props, open(os.path.join(OUT, 'props.json'), 'w'), default=str); props = json.load(open(os.path.join(OUT, 'props.json')))
adj = analyze.adjacency(props, n)
def draw(p, shots): return (rng.random((shots, len(p))) < p).astype(np.int8)
def spec_bits(c, phi, st, shots):   # coherence c, phase phi (rad) -> bits for X or Y reading
    e = c * (np.cos(phi) if st == 'XX' else np.sin(phi)); return draw((1 - e) / 2, shots)
PLANT = {}
def make(run, chip='ibm_fez'):
    pubs, meta, shots, info = runs.BUILD[run](b, chip)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; z = {}
    for i, m in enumerate(meta):
        if run == 'distance':
            E, S = info['E_%s' % m['layout']], info['S_%s' % m['layout']]; dd = analyze.bfs(adj, E); d = np.array([dd.get(q, 99) for q in S])
            kick = -np.radians(3.0) * np.exp(-(d - 1) / 2.0) if m['layout'] == 'ball' else -np.radians(2.0) * np.ones(len(S))
            if NULL: kick = 0 * kick
            phi = (kick * m['N']) if m['mode'] == 'M' else 0 * d; z['fin_%d' % i] = spec_bits(0.8, phi, m['setting'], shots)
        if run == 'count':
            N = m['N']; cD = 0.9 * np.exp(-N * 0.002); c = cD * (0.995 ** N if m['mode'] == 'M' and not NULL else 1); phi = np.radians(-2.0) * N if m['mode'] == 'M' and not NULL else 0.0
            z['fin_%d' % i] = spec_bits(c, phi * np.ones(len(info['S'])), m['setting'], shots)
        if run == 'order':
            G = info['groups']; p = np.zeros(n) + 0.01
            if m['kind'] == 'alone' and not NULL: p += 0.002 * m['slot']
            if m['kind'] == 'order' and not NULL:
                for k, g in enumerate(m['order']): p[G[g]] += 0.002 * k + (0, 0.0015, 0.003)[k]
            z['c_%d' % i] = draw(p, shots)
        if run == 'gates':
            p = np.zeros(n) + 0.02
            if 'ngates' in m: p = p + (0 if NULL else 0.0005 * (m['reps'] - 1) * m['ngates'])
            else: p = np.zeros(n) + np.exp(-m['ref_wait_us'] / 150)
            z['c_%d' % i] = draw(p, shots)
        if run == 'neighbour':
            s = 1e-4 if NULL else {'0': 1e-4, '1': 1.5e-4, '+': 1.25e-4}[m['arm']]; z['t_%d' % i] = draw(np.zeros(len(info['targets'])) + 0.005 + s * m['m'], shots)
        if run == 'revival':
            t = m['tau']; c = 0.9 * np.exp(-t / 40) * (abs(np.cos(np.pi * t / 68)) if m['arm'] == '+' and not NULL else 1) + 0.02
            st = 'XX' if m['setting'] == 'X' else 'YY'; z['a_%d' % i] = spec_bits(c, np.zeros(len(info['A'])) + 0.3, st, shots)
    rec = dict(run=run, chip=chip, n_qubits=n, job_id='SYNTH', meta=meta, info=info, shots=shots)
    json.dump(rec, open(os.path.join(OUT, 'job_%s_%s.json' % (run, chip)), 'w'), default=str)
    json.dump(props, open(os.path.join(OUT, 'props_%s_%s.json' % (run, chip)), 'w'))
    np.savez_compressed(os.path.join(OUT, 'raw_%s_%s.npz' % (run, chip)), **z)
for r in (os.environ.get('RUNS') or ','.join(runs.BUILD)).split(','): make(r)
if NULL: print('NULL: no effect planted anywhere, seed %d' % SEEDN)
else: print('PLANTED: distance ball kick -3.0 x exp(-(d-1)/2) deg (falls off: R1.3 BELOW, R1.4 ABOVE), spread -2.0 flat; count kick -2.0 deg flat (R2.2, R2.3 CONSISTENT),'
      ' per-readout product loss 0.005 constant (R2.7 CONSISTENT, R2.6 ~0.005); order slot excess 0 / 0.0015 / 0.003, waiting +0.004 slot3-slot1, together 0;'
      ' gates 0.0005 per extra-gate-step = 0.00025 per extra gate; neighbour 1e-4 / 1.5e-4 / 1.25e-4 (R5.4 CONSISTENT); revival: beat P=68 in +, plain decay in 0 and 1')
analyze.run_all(OUT)
