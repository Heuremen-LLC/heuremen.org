"""Reach run, 22 Oct 2026 (Heuremen LLC, Wayfinder and Bones). PLAN_REACH.md, sha256 in PLAN_REACH.sha256. Where does the readout kick stop?
Dry-run by default: builds, picks the layouts from the live coupling map, prints the untouched-qubit counts at every exact distance, transpiles
against the live chip, estimates seconds, submits nothing. --fire submits (Darkwave Disco instance). Usage: python reach.py CHIP [--fire] [--shots S]
Design: K = 3 packed layouts in one job. Each layout is a breadth-first ball of 39 measured qubits around a centre; the centres are chosen from the
coupling map BEFORE any data by an exhaustive search over every triple of centres (rule in choose_layouts). Per layout: idle (D) vs readout (M),
N = 8 readouts, untouched qubits on the equator read in X and Y. One shared N = 0 baseline (all 156 qubits on the equator, read in X and Y).
Every circuit twice per job (copy 0, copy 1), shuffled together."""
import os, sys, json, time, random, itertools, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import transpile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
sys.path.append(os.path.join(ROOT, 'control_run')); sys.path.append(os.path.join(ROOT, 'deciding_runs'))
import run as base
from runs import edges, adjacency, bfs, set_L
SEED = 20261022; K_LAYOUTS = 3; N_MEAS = 39; N_READ = 8; SHOTS = 4000
BINS = (3, 4, 5, 6, 7, 8)          # exact distance to the nearest measured qubit; 8 means 8 or more
def binof(d): return min(d, 8)
def ball(adj, c):
    order = sorted(bfs(adj, [c]).items(), key=lambda kv: (kv[1], kv[0])); return sorted(q for q, _ in order[:N_MEAS])
def shells(adj, n, E):
    dd = bfs(adj, E); Es = set(E); S = [q for q in range(n) if q not in Es]; return S, [int(dd.get(q, 99)) for q in S]
def choose_layouts(b):
    """Every qubit is a candidate centre. For every set of K centres, count per bin (3, 4, 5, 6, 7, 8+) the DISTINCT physical untouched qubits that sit
    at that exact distance in at least one layout, and the pooled (layout, qubit) units. Pick the set that maximises, in order: the smallest distinct
    count over bins 3-7; the smallest distinct count over all six bins; the smallest pooled count over bins 3-7; the total pooled count over the six
    bins; then a seeded random rank (SEED) as the last tie-break. Coupling map only; no data."""
    n = b.num_qubits; adj = adjacency(b); bits = {}
    for c in range(n):
        S, d = shells(adj, n, ball(adj, c)); m = [0] * len(BINS)
        for q, dq in zip(S, d):
            if dq >= 3: m[BINS.index(binof(dq))] |= 1 << q
        bits[c] = m
    rank = list(range(n)); random.Random(SEED).shuffle(rank); rank = {c: i for i, c in enumerate(rank)}
    pc = {c: [bin(x).count('1') for x in bits[c]] for c in range(n)}
    best, bkey = None, None
    for cs in itertools.combinations(range(n), K_LAYOUTS):
        u = [bin(bits[cs[0]][k] | bits[cs[1]][k] | bits[cs[2]][k]).count('1') for k in range(len(BINS))]
        p = [pc[cs[0]][k] + pc[cs[1]][k] + pc[cs[2]][k] for k in range(len(BINS))]
        key = (min(u[:5]), min(u), min(p[:5]), sum(p), -sum(rank[c] for c in cs))
        if bkey is None or key > bkey: best, bkey = cs, key
    return list(best)
def build(b, chip, shots=SHOTS, verbose=True):
    n = b.num_qubits; adj = adjacency(b); centres = choose_layouts(b); lays = []
    for c in centres:
        E = ball(adj, c); S, d = shells(adj, n, E); lays.append(dict(centre=c, E=E, S=S, d=d))
    L = set_L(b, lays[0]['E'][0]); Ls = {int(round(b.target['measure'][(q,)].duration / b.target.dt)) for q in range(n)}
    if verbose:
        print('%s: centres %s; readout length L = %d dt (all qubits: %s)' % (chip, centres, L, sorted(Ls)))
        hdr = '  untouched qubits at exact distance  ' + ' '.join('%5s' % ('d=%d' % k if k < 8 else 'd8+') for k in (1, 2) + BINS)
        print(hdr)
        for lay in lays:
            d = np.array(lay['d']); print('  layout centre %-4d                  ' % lay['centre'] + ' '.join('%5d' % (np.sum(d == k) if k < 8 else np.sum(d >= 8)) for k in (1, 2) + BINS))
        pooled = [sum(int(np.sum(np.array(l['d']) == k) if k < 8 else np.sum(np.array(l['d']) >= 8)) for l in lays) for k in (1, 2) + BINS]
        distinct = [len({q for l in lays for q, dq in zip(l['S'], l['d']) if (dq == k if k < 8 else dq >= 8)}) for k in (1, 2) + BINS]
        print('  pooled (layout, qubit) units        ' + ' '.join('%5d' % x for x in pooled))
        print('  distinct physical qubits            ' + ' '.join('%5d' % x for x in distinct))
    pubs, meta = [], []
    for st in ('XX', 'YY'):
        pubs.append(base.spectator_circuit(n, [], [], list(range(n)), [], 'D', 0, st, '0', False)); meta.append(dict(layout=-1, mode='D', N=0, setting=st))
    for k, lay in enumerate(lays):
        for mode in ('D', 'M'):
            for st in ('XX', 'YY'):
                pubs.append(base.spectator_circuit(n, lay['E'], [], lay['S'], lay['E'], mode, N_READ, st, '0', False)); meta.append(dict(layout=k, mode=mode, N=N_READ, setting=st))
    info = dict(layouts=lays, L_dt=L, bins=list(BINS), n_meas=N_MEAS, N=N_READ, coupling_edges=[list(e) for e in edges(b)])
    return pubs, meta, shots, info
if __name__ == '__main__':
    CHIP = sys.argv[1]; fire = '--fire' in sys.argv; shots = int(sys.argv[sys.argv.index('--shots') + 1]) if '--shots' in sys.argv else SHOTS
    s = base.service(); rem = s.usage().get('usage_remaining_seconds') or 0; b = s.backend(CHIP)
    pubs, meta, shots, info = build(b, CHIP, shots)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED + 7).shuffle(order); pubs = [pubs[i] for i in order]; meta = [meta[i] for i in order]
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(shots * np.sum(dur + 290e-6))
    print('reach %s: %d pubs, %d shots each, circuits %.1f-%.1f us, total shots %d, ESTIMATED %.0f s, remaining %s, queue %s'
          % (CHIP, len(tq), shots, dur.min() * 1e6, dur.max() * 1e6, shots * len(tq), est, rem, b.status().pending_jobs), flush=True)
    if not fire: sys.exit()
    if rem < est + 3: print('STOP: estimated %.0f s, %.0f s left' % (est, rem)); sys.exit()
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(HERE + '/props_reach_%s.json' % CHIP, 'w'), default=str)
    from qiskit_ibm_runtime import SamplerV2
    job = SamplerV2(mode=b).run(tq, shots=shots)
    rec = dict(run='reach', chip=CHIP, n_qubits=b.num_qubits, job_id=job.job_id(), submitted=time.strftime('%Y-%m-%dT%H:%M:%S%z'), meta=meta, order=order,
               shots=shots, info=info, plan_sha256=open(HERE + '/PLAN_REACH.sha256').read().split()[0], code_sha256=open(HERE + '/CODE_REACH.sha256').read())
    json.dump(rec, open(HERE + '/job_reach_%s.json' % CHIP, 'w'), default=str); print('SUBMITTED reach', CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        for name in pr.data: out['%s_%d' % (name, i)] = getattr(pr.data, name).array
    np.savez_compressed(HERE + '/raw_reach_%s.npz' % CHIP, **out)
    print('DONE reach %s: remaining now %s' % (CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
