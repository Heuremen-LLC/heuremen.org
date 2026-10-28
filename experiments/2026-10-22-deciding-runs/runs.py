"""Deciding runs, 22 Oct 2026 (Heuremen LLC). PLAN.md, sha256 in PLAN.sha256. Dry-run by default: builds, transpiles against the live chip,
estimates seconds, submits nothing. --fire submits. Darkwave Disco instance. Usage: python runs.py RUN CHIP [--fire]
RUN in: distance, count, order, gates, neighbour, revival. Every circuit twice per job (copy 0, copy 1), shuffled together."""
import os, sys, json, time, random, warnings; warnings.filterwarnings('ignore')
from collections import deque
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import SamplerV2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.append(os.path.join(HERE, '..', 'control_run'))
import run as base
SEED = 20261022
def edges(b): return sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()})
def adjacency(b):
    adj = {q: set() for q in range(b.num_qubits)}
    for a, c in edges(b): adj[a].add(c); adj[c].add(a)
    return adj
def bfs(adj, src):
    d = {s: 0 for s in src}; q = deque(src)
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in d: d[v] = d[u] + 1; q.append(v)
    return d
def set_L(b, q): base.L_DT = int(round(b.target['measure'][(q,)].duration / b.target.dt)); return base.L_DT
def rot(qc, q, basis):
    if basis == 'X': qc.h(q)
    elif basis == 'Y': qc.sdg(q); qc.h(q)
# ---------------- run 1: distance. 39 measured qubits packed in one region (breadth-first ball from qubit 0), the other 117 untouched at
# distances 1 to the far side; plus the spread layout of the control run in the same job (like for like).
def build_distance(b, chip):
    n = b.num_qubits; adj = adjacency(b); order = sorted(bfs(adj, [0]).items(), key=lambda kv: (kv[1], kv[0]))
    Eball = sorted([q for q, _ in order[:39]]); Sball = [q for q in range(n) if q not in set(Eball)]
    Espread = sorted(json.load(open(base.ROLES % chip))['emitters']); Sspread = [q for q in range(n) if q not in set(Espread)]
    set_L(b, Eball[0]); pubs, meta = [], []
    for lay, E, S in (('ball', Eball, Sball), ('spread', Espread, Sspread)):
        for st in ('XX', 'YY'): pubs.append(base.spectator_circuit(n, E, [], S, E, 'D', 0, st, '0', False)); meta.append(dict(layout=lay, mode='D', N=0, setting=st))
        for mode in ('D', 'M'):
            for st in ('XX', 'YY'): pubs.append(base.spectator_circuit(n, E, [], S, E, mode, 8, st, '0', False)); meta.append(dict(layout=lay, mode=mode, N=8, setting=st))
    return pubs, meta, 4000, dict(E_ball=Eball, S_ball=Sball, E_spread=Espread, S_spread=Sspread, L_dt=base.L_DT)
# ---------------- run 2: readout count. Control-run layout, idle vs readout at N = 1, 3, 9, 27, 54, 81.
COUNT_N = (1, 3, 9, 27, 54, 81)
def build_count(b, chip):
    n = b.num_qubits; E = sorted(json.load(open(base.ROLES % chip))['emitters']); S = [q for q in range(n) if q not in set(E)]
    set_L(b, E[0]); pubs, meta = [], []
    for st in ('XX', 'YY'): pubs.append(base.spectator_circuit(n, E, [], S, E, 'D', 0, st, '0', False)); meta.append(dict(mode='D', N=0, setting=st))
    for N in COUNT_N:
        for mode in ('D', 'M'):
            for st in ('XX', 'YY'): pubs.append(base.spectator_circuit(n, E, [], S, E, mode, N, st, '0', False)); meta.append(dict(mode=mode, N=N, setting=st))
    return pubs, meta, 1000, dict(E=E, S=S, L_dt=base.L_DT)
# ---------------- run 3: read order. All qubits in 0, split into groups A, B, C (seeded). Three readout slots of length L. Every order of the
# three groups; each group read ALONE at each slot (everyone else idle, nobody else read); all three read together at slot 1.
def build_order(b, chip):
    n = b.num_qubits; Q = list(range(n)); random.Random(SEED).shuffle(Q); G = {g: sorted(Q[i::3]) for i, g in enumerate('ABC')}
    L = set_L(b, 0)
    def circ(slots):   # slots: list of 3 lists of group names read in that slot
        qc = QuantumCircuit(n, n)
        for s in slots:
            rd = [q for g in s for q in G[g]]; idle = [q for q in range(n) if q not in set(rd)]
            if rd: qc.measure(rd, rd)
            if idle: qc.delay(L, idle, unit='dt')
            qc.barrier()
        return qc
    pubs, meta = [], []
    import itertools
    for perm in itertools.permutations('ABC'): pubs.append(circ([[perm[0]], [perm[1]], [perm[2]]])); meta.append(dict(kind='order', order=''.join(perm)))
    for g in 'ABC':
        for k in range(3):
            sl = [[], [], []]; sl[k] = [g]; pubs.append(circ(sl)); meta.append(dict(kind='alone', group=g, slot=k))
    pubs.append(circ([['A', 'B', 'C'], [], []])); meta.append(dict(kind='together'))
    return pubs, meta, 3000, dict(groups=G, L_dt=L)
def gate_circuit(n, ngates, reps, ref_wait=None):
    """The control run's gate_circuit with a barrier after EVERY X, so X X X is not compiled down to X (found in the dry run)."""
    qc = QuantumCircuit(n, n); Q = list(range(n))
    if ref_wait is not None:
        qc.x(Q); qc.barrier(); qc.delay(ref_wait, Q, unit='us'); qc.barrier(); qc.measure(Q, Q); return qc
    for _ in range(ngates):
        qc.delay(base.WAIT_US, Q, unit='us'); qc.barrier()
        for _ in range(reps): qc.x(Q); qc.barrier()
    qc.measure(Q, Q); return qc
# ---------------- run 4: gates. The control run's gate arm unchanged: X trains of 2, 6, 10 (30 us apart) as X or as X X X on identical timing;
# references X then wait 30, 90, 300 us.
def build_gates(b, chip):
    base.ARMS = ['G']; base.CHIP = chip; base.gate_circuit = gate_circuit
    pubs, meta, *_ = base.build(b, gate_only=True)
    keep = [i for i, m in enumerate(meta) if m['copy'] == 0]; return [pubs[i] for i in keep], [{k: v for k, v in meta[i].items() if k != 'copy'} for i in keep], 1000, {}
# ---------------- run 5: neighbour state and gate error. Targets at least 3 apart; their neighbours prepared in 0, 1 or + (never read);
# each target runs m = 0, 40, 80, 160 X gates back to back and is read in Z.
NB_M = (0, 40, 80, 160)
def build_neighbour(b, chip):
    n = b.num_qubits; adj = adjacency(b); T = []
    for q in range(n):
        if all(bfs(adj, [q]).get(t, 99) >= 3 for t in T): T.append(q)
    Nb = sorted({v for t in T for v in adj[t]})
    pubs, meta = [], []
    for arm in ('0', '1', '+'):
        for m in NB_M:
            qc = QuantumCircuit(n); qc.add_register(ClassicalRegister(len(T), 't'))
            if arm == '1': qc.x(Nb)
            if arm == '+': qc.h(Nb)
            qc.barrier()
            for _ in range(m): qc.x(T); qc.barrier(T)
            qc.barrier(); qc.measure(T, range(len(T))); pubs.append(qc); meta.append(dict(arm=arm, m=m))
    return pubs, meta, 2000, dict(targets=T, neighbours=Nb)
# ---------------- run 6: revival. Two-colour the coupling graph: colour A runs a free Ramsey clock (SX, wait tau, read X and Y); colour B is
# prepared in 0, 1 or + and never read. tau = 2 to 150 us in steps of 4.
TAUS = tuple(range(2, 151, 4))
def build_revival(b, chip):
    n = b.num_qubits; adj = adjacency(b); col = {}
    for s in range(n):
        if s in col: continue
        col[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in col: col[v] = 1 - col[u]; q.append(v)
                elif col[v] == col[u]: raise SystemExit('coupling graph not two-colourable on %s' % chip)
    A = [q for q in range(n) if col[q] == 0]; B = [q for q in range(n) if col[q] == 1]
    pubs, meta = [], []
    for arm in ('0', '1', '+'):
        for tau in TAUS:
            for st in 'XY':
                qc = QuantumCircuit(n); qc.add_register(ClassicalRegister(len(A), 'a'))
                if arm == '1': qc.x(B)
                if arm == '+': qc.h(B)
                qc.sx(A); qc.barrier(); qc.delay(tau, A + B, unit='us'); qc.barrier()
                for q in A: rot(qc, q, st)
                qc.measure(A, range(len(A))); pubs.append(qc); meta.append(dict(arm=arm, tau=tau, setting=st))
    return pubs, meta, 500, dict(A=A, B=B)
BUILD = dict(distance=build_distance, count=build_count, order=build_order, gates=build_gates, neighbour=build_neighbour, revival=build_revival)
if __name__ == '__main__':
    RUN, CHIP = sys.argv[1], sys.argv[2]; fire = '--fire' in sys.argv
    s = base.service(); rem = s.usage().get('usage_remaining_seconds') or 0; b = s.backend(CHIP)
    pubs, meta, shots, info = BUILD[RUN](b, CHIP)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED + 7).shuffle(order); pubs = [pubs[i] for i in order]; meta = [meta[i] for i in order]
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(shots * np.sum(dur + 290e-6))
    print('%s %s: %d pubs, %d shots each, circuits %.1f-%.1f us, total shots %d, ESTIMATED %.0f s, remaining %s, queue %s'
          % (RUN, CHIP, len(tq), shots, dur.min() * 1e6, dur.max() * 1e6, shots * len(tq), est, rem, b.status().pending_jobs), flush=True)
    if not fire: sys.exit()
    if rem < est + 3: print('STOP: estimated %.0f s, %.0f s left' % (est, rem)); sys.exit()
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(HERE + '/props_%s_%s.json' % (RUN, CHIP), 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=shots)
    rec = dict(run=RUN, chip=CHIP, n_qubits=b.num_qubits, job_id=job.job_id(), submitted=time.strftime('%Y-%m-%dT%H:%M:%S%z'), meta=meta, order=order, shots=shots,
               info=info, plan_sha256=open(HERE + '/PLAN.sha256').read().split()[0])
    json.dump(rec, open(HERE + '/job_%s_%s.json' % (RUN, CHIP), 'w'), default=str); print('SUBMITTED', RUN, CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        for name in pr.data: out['%s_%d' % (name, i)] = getattr(pr.data, name).array
    np.savez_compressed(HERE + '/raw_%s_%s.npz' % (RUN, CHIP), **out)
    print('DONE %s %s: remaining now %s' % (RUN, CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
