"""Control run (PREREG.md, sha256 in PREREG.sha256): one job per chip that tests the ordinary causes of every readout claim.
Dry-run by default (builds, transpiles, estimates seconds; nothing submitted). --fire submits. Darkwave Disco instance."""
import os, sys, json, time, random, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
H = os.path.expanduser('~/bones-forge'); HERE = os.path.dirname(os.path.abspath(__file__))
ROLES = H + '/heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_%s.json'
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
CHIP = arg('--chip', 'ibm_kingston'); SHOTS = int(arg('--shots', 1000)); SEED = 20261021
SETTINGS = ['ZZ', 'XX', 'YY', 'XY', 'YX']
WAIT_US = 30
def service():
    tok = open(H + '/.ibm_quantum_key_3').read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == 'Darkwave Disco'][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)
def pairs_for(b, S):
    Sset = set(S); used = set(); P = []
    for a, c in sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}):
        if a in Sset and c in Sset and a not in used and c not in used: P.append((a, c)); used |= {a, c}
    return P
def rot(qc, q, basis):
    if basis == 'X': qc.h(q)
    elif basis == 'Y': qc.sdg(q); qc.h(q)
def spectator_circuit(n, E, P, SG, K_set, mode, N, st, e_init, e_final_read, leak_flip=False):
    """E = all 39 emitters. K_set = the emitters actually measured in M modes (the rest idle exactly L).
    mode 'D' idles K_set for L per step, 'M' measures K_set each step. e_init '0' or '1'."""
    regs = [ClassicalRegister(len(K_set), 'm%d' % k) for k in range(N)]
    fin = ClassicalRegister(2 * len(P) + len(SG), 'fin'); qc = QuantumCircuit(n)
    for r in regs: qc.add_register(r)
    qc.add_register(fin)
    efin = None
    if e_final_read: efin = ClassicalRegister(len(E), 'efin'); qc.add_register(efin)
    if e_init == '1': qc.x(E)
    for a, c in P: qc.h(a); qc.cx(a, c)
    for q in SG: qc.sx(q)
    qc.barrier()
    for k in range(N):
        if mode == 'M':
            qc.measure(K_set, regs[k])
            for q in E:
                if q not in K_set: qc.delay(L_DT, q, unit='dt')
        else:
            for q in E: qc.delay(L_DT, q, unit='dt')
        qc.barrier()
    for a, c in P: rot(qc, a, st[0]); rot(qc, c, st[1])
    for q in SG: rot(qc, q, st[0])
    for i, (a, c) in enumerate(P): qc.measure(a, fin[2 * i]); qc.measure(c, fin[2 * i + 1])
    for i, q in enumerate(SG): qc.measure(q, fin[2 * len(P) + i])
    if e_final_read:
        if leak_flip: qc.x(E)
        qc.measure(E, efin)
    return qc
def gate_circuit(n, ngates, reps, ref_wait=None):
    """All qubits from |0>. Train: ngates x [wait 30 us, X^reps]; read right after the last gate.
    reps=1 is the original design; reps=3 applies X X X (same logical X, three times the gates) on identical timing.
    Reference: X, then wait ref_wait us, read."""
    qc = QuantumCircuit(n, n); Q = list(range(n))
    if ref_wait is not None:
        qc.x(Q); qc.barrier(); qc.delay(ref_wait, Q, unit='us'); qc.barrier(); qc.measure(Q, Q); return qc
    for _ in range(ngates):
        qc.delay(WAIT_US, Q, unit='us'); qc.barrier()
        for _ in range(reps): qc.x(Q)
        qc.barrier()
    qc.measure(Q, Q); return qc
# shots per copy, per arm (two copies each). Power set from the 21 Oct data:
#  S/C singles: 117 qubits, coherence width ~0.003 at 2000 total shots (effect 0.005-0.013); kick width ~0.15 deg (effect ~2 deg)
#  P pairs: 48 pairs, concurrence cost ~0.005-0.015 needs ~8000 total shots to resolve on fez; N kept at 2 and 4 (marrakesh saturates by 8)
#  G gate: 800 shots gave slope width ~0.0003; X vs XXX must separate -0.0014 from -0.0042
ARMS = arg('--arms', 'S,C,P,L').split(',')   # the wild claims and their controls; G (gate) and B (measured qubit) on request
SHOTS_ARM = {'S': 1000, 'C': 1000, 'P': 4000, 'B': 1000, 'L': 1000, 'G': 1000}
def build(b, gate_only=False):
    global L_DT
    n = b.num_qubits; pubs, meta = [], []
    def add(qc, **m): pubs.append(qc); meta.append(m)
    if not gate_only:
        r = json.load(open(ROLES % CHIP)); E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
        P = pairs_for(b, S); paired = {q for pr in P for q in pr}; SG = [q for q in S if q not in paired]
        L_DT = int(round(b.target['measure'][(E[0],)].duration / b.target.dt))
        rng = random.Random(SEED); K10 = sorted(rng.sample(E, 10)); K1 = [K10[0]]
        # S: all 117 untouched qubits unpaired on the equator; 39 measured; idle N=0 and idle/readout at N=2,4,8; read X and Y
        for st in ('XX', 'YY'): add(spectator_circuit(n, E, [], S, E, 'D', 0, st, '0', False), arm='S', mode='D', N=0, K=39, setting=st)
        for N in (2, 4, 8):
            for mode in ('D', 'M'):
                for st in ('XX', 'YY'): add(spectator_circuit(n, E, [], S, E, mode, N, st, '0', False), arm='S', mode=mode, N=N, K=39, setting=st)
        # C: dose, 1 and 10 of the 39 measured (rest idle exactly L), N=8, same unpaired spectators
        for name, K in (('K1', K1), ('K10', K10)):
            for mode in ('D', 'M'):
                for st in ('XX', 'YY'): add(spectator_circuit(n, E, [], S, K, mode, 8, st, '0', False), arm='C', mode=mode, N=8, K=len(K), setting=st)
        # P: 48 Bell pairs (+21 singles), idle N=0 and idle/readout at N=2,4, five tomography settings
        for st in SETTINGS: add(spectator_circuit(n, E, P, SG, E, 'D', 0, st, '0', False), arm='P', mode='D', N=0, K=39, setting=st)
        for N in (2, 4):
            for mode in ('D', 'M'):
                for st in SETTINGS: add(spectator_circuit(n, E, P, SG, E, mode, N, st, '0', False), arm='P', mode=mode, N=N, K=39, setting=st)
        # B: measured qubits prepared in 1, idle vs readout at N=4,8, final readout of every measured qubit (direct, no decay model)
        for N in (4, 8):
            for mode in ('D', 'M'): add(spectator_circuit(n, E, [], S, E, mode, N, 'XX', '1', True), arm='B', mode=mode, N=N, K=39, setting='XX')
        # L: leakage. E in 0, 8 readouts, X on E, final readout; reference: E in 1, 1 readout, X, final readout
        add(spectator_circuit(n, E, [], S, E, 'M', 8, 'XX', '0', True, leak_flip=True), arm='L', mode='M', N=8, K=39, setting='XX', e_init='0')
        add(spectator_circuit(n, E, [], S, E, 'M', 1, 'XX', '1', True, leak_flip=True), arm='L', mode='M', N=1, K=39, setting='XX', e_init='1')
    else:
        E = S = P = SG = K1 = K10 = []; L_DT = None
    # G: gate decider on every qubit: n = 2, 6, 10 trains of X vs XXX on identical timing; references X then wait 30/90/300 us
    for ng in (2, 6, 10):
        for reps in (1, 3): add(gate_circuit(n, ng, reps), arm='G', ngates=ng, reps=reps)
    for w in (30, 90, 300): add(gate_circuit(n, 0, 1, ref_wait=w), arm='G', ref_wait_us=w)
    keep = [i for i, m in enumerate(meta) if m['arm'] in ARMS]; pubs = [pubs[i] for i in keep]; meta = [meta[i] for i in keep]
    # every pub twice in the same job (copy 0 and copy 1), shuffled together: the copy-to-copy gap is measured, not assumed zero
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED + 1).shuffle(order)
    return [pubs[i] for i in order], [meta[i] for i in order], order, E, S, P, SG, K1, K10

if __name__ == '__main__':
    fire = '--fire' in sys.argv; s = service(); rem = s.usage().get('usage_remaining_seconds') or 0
    print('remaining seconds:', rem)
    GATE_ONLY = '--gate-only' in sys.argv
    b = s.backend(CHIP); pubs, meta, order, E, S, P, SG, K1, K10 = build(b, GATE_ONLY)
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq])
    shots = np.array([SHOTS_ARM[m['arm']] for m in meta]); est = float(np.sum(shots * (dur + 290e-6)))   # per-shot model from two 21 Oct jobs only (70000 shots ~ 22-24 s; 36000 ~ 12 s); circuits here run to 300 us, past the 90 us they were fitted on
    arms = {}
    for m in meta: arms[m['arm']] = arms.get(m['arm'], 0) + 1
    print('%s: %d pubs %s, E %d, S %d, pairs %d, singles %d, L=%s dt, circuits %.1f-%.1f us, total shots %d, ESTIMATED %.0f s, queue %s'
          % (CHIP, len(tq), arms, len(E), len(S), len(P), len(SG), L_DT, dur.min() * 1e6, dur.max() * 1e6, int(shots.sum()), est, b.status().pending_jobs))
    if not fire: sys.exit()
    if rem < est + 3: print('STOP: estimated %.0f s, %.0f s left' % (est, rem)); sys.exit()
    props = b.properties()   # IBM's stated gate/readout errors and T1/T2 at the moment of firing, saved with the job
    json.dump(props.to_dict() if props is not None else None, open(HERE + '/backend_properties_%s.json' % CHIP, 'w'), default=str)
    job = SamplerV2(mode=b).run([(c, None, int(k)) for c, k in zip(tq, shots)])
    rec = {'chip': CHIP, 'job_id': job.job_id(), 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'meta': meta, 'order': order, 'E': E, 'S': S,
           'pairs': P, 'singles': SG, 'K1': K1, 'K10': K10, 'L_dt': L_DT, 'shots': [int(k) for k in shots], 'gate_only': GATE_ONLY, 'arms': ARMS, 'prereg_sha256': open(HERE + '/PREREG.sha256').read().split()[0]}
    json.dump(rec, open(HERE + '/job_%s.json' % CHIP, 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        for name in pr.data:
            out['%s_%d' % (name, i)] = getattr(pr.data, name).array
    np.savez_compressed(HERE + '/raw_%s.npz' % CHIP, **out)
    print('DONE %s: remaining now %s' % (CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
