"""Flat-beat prediction run on fresh phoenix pairs, 22 Oct 2026 (built from zz.py of the one-neighbour run) (Heuremen LLC, Wayfinder and Bones). PLAN_FLAT.md, sha256 in PLAN_FLAT.sha256.
Dry run by default (builds, transpiles against the live chip, estimates seconds, submits nothing). --fire submits.
Usage: python zz.py CHIP [--fire]
Layout: triples (B1, A, B2), B1 and B2 coupled to A. Every other neighbour of A is idle in 0 (no A of another triple and no B of
another triple is coupled to A), so A can beat only against B1 and B2. A runs a free clock (SX, wait tau, read X or Y); B1, B2 are
read in Z at the end. Arms: Z (B1, B2 in 0), P1 (B1 in +), P2 (B2 in +), P12 (both in +), O1 (B1 in 1), EZ and EP1 (as Z and P1 with a
Hahn echo X on A at tau/2). Every circuit twice per job (copy 0, copy 1), shuffled together."""
import os, sys, json, time, random, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import SamplerV2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.append(os.path.join(HERE, '..', 'control_run'))
import run as base
SEED = 20261024; SHOTS = 2000
ARMS = ('Z', 'P1', 'P2')
TAU_US = {'ibm_phoenix': [2 + 4 * k for k in range(41)]}   # 2-162 us
def adjacency(b):
    adj = {q: set() for q in range(b.num_qubits)}
    for a, c in b.coupling_map.get_edges(): adj[a].add(c); adj[c].add(a)
    return adj
def qprops(props, q):
    d = {p['name']: p['value'] for p in props['qubits'][q]}; return d.get('readout_error', 1.0), d.get('T1', 0.0), d.get('T2', 0.0)
def layout(b, props):
    """Greedy, seeded. A needs readout error < 0.05 and T2 > 40 us; B1, B2 readout error < 0.05."""
    adj = adjacency(b); n = b.num_qubits; rng = random.Random(SEED); order = list(range(n)); rng.shuffle(order)
    ok_a = lambda q: qprops(props, q)[0] < 0.05 and qprops(props, q)[2] > 40
    ok_b = lambda q: qprops(props, q)[0] < 0.05
    prev = json.load(open(os.path.join(HERE, '..', 'zz_revival', 'job_zz_ibm_phoenix.json')))['info']['triples']
    old_A = {t[1] for t in prev}   # fresh pairs: no A of the 22 Oct one-neighbour run is an A here
    used, As, trip = set(), set(), []
    for a in order:
        if a in used or a in old_A or not ok_a(a) or adj[a] & used: continue
        cand = sorted(q for q in adj[a] if q not in used and ok_b(q) and not (adj[q] & As)); rng.shuffle(cand)
        if len(cand) < 2: continue
        b1, b2 = cand[0], cand[1]
        trip.append((b1, a, b2)); used |= {a, b1, b2}; As.add(a)
    return trip
def tau_dt(b, us):
    g = 16; dt = b.target.dt; return int(round(us * 1e-6 / dt / g)) * g
def rot(qc, q, basis):
    if basis == 'X': qc.h(q)
    elif basis == 'Y': qc.sdg(q); qc.h(q)
def build(b, chip, props):
    trip = layout(b, props); n = b.num_qubits
    A = [t[1] for t in trip]; B1 = [t[0] for t in trip]; B2 = [t[2] for t in trip]
    pubs, meta = [], []
    for arm in ARMS:
        for tus in TAU_US[chip]:
            T = tau_dt(b, tus)
            for st in 'XY':
                qc = QuantumCircuit(n); ra = ClassicalRegister(len(A), 'a'); r1 = ClassicalRegister(len(A), 'b1'); r2 = ClassicalRegister(len(A), 'b2')
                qc.add_register(ra); qc.add_register(r1); qc.add_register(r2)
                if arm in ('P1', 'P12', 'EP1'): qc.h(B1)
                if arm in ('P2', 'P12'): qc.h(B2)
                if arm == 'O1': qc.x(B1)
                qc.sx(A); qc.barrier()
                if arm in ('EZ', 'EP1'):
                    h = (T // 2 // 16) * 16
                    qc.delay(h, A + B1 + B2, unit='dt'); qc.barrier(); qc.x(A); qc.barrier(); qc.delay(T - h, A + B1 + B2, unit='dt')
                else:
                    qc.delay(T, A + B1 + B2, unit='dt')
                qc.barrier()
                for q in A: rot(qc, q, st)
                qc.measure(A, ra); qc.measure(B1, r1); qc.measure(B2, r2)
                pubs.append(qc); meta.append(dict(arm=arm, tau_us=tus, tau_dt=T, setting=st))
    return pubs, meta, dict(triples=trip)
if __name__ == '__main__':
    CHIP = sys.argv[1]; fire = '--fire' in sys.argv
    s = base.service(); rem = s.usage().get('usage_remaining_seconds') or 0; b = s.backend(CHIP)
    props = b.properties().to_dict()
    pubs, meta, info = build(b, CHIP, props)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED + 7).shuffle(order); pubs = [pubs[i] for i in order]; meta = [meta[i] for i in order]
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    nx = sum(1 for c in tq[:50] for ins in c.data if ins.operation.name == 'x')
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6))
    print('zz %s: %d triples, %d pubs, %d shots each, circuits %.1f-%.1f us, ESTIMATED %.0f s, remaining %s, queue %s, x ops in first 50: %d'
          % (CHIP, len(info['triples']), len(tq), SHOTS, dur.min() * 1e6, dur.max() * 1e6, est, rem, b.status().pending_jobs, nx), flush=True)
    if not fire: sys.exit()
    if rem < est + 3: print('STOP: estimated %.0f s, %.0f s left' % (est, rem)); sys.exit()
    json.dump(props, open(HERE + '/props_flat_%s.json' % CHIP, 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = dict(run='flat', chip=CHIP, n_qubits=b.num_qubits, job_id=job.job_id(), submitted=time.strftime('%Y-%m-%dT%H:%M:%S%z'), meta=meta, order=order,
               shots=SHOTS, info=info, edges=sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}), plan_sha256=open(HERE + '/PLAN_FLAT.sha256').read().split()[0])
    json.dump(rec, open(HERE + '/job_flat_%s.json' % CHIP, 'w'), default=str); print('SUBMITTED zz', CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        for name in pr.data: out['%s_%d' % (name, i)] = getattr(pr.data, name).array
    np.savez_compressed(HERE + '/raw_flat_%s.npz' % CHIP, **out)
    print('DONE zz %s: remaining now %s' % (CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
