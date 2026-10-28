"""run_len.py -- two readout lengths on one chip (PLAN_LENGTH.md). Usage: run_len.py CHIP [--fire] [--shots N]
Same 39 measured / 117 untouched roles as the control run (21 Oct). Every readout step is the chip's `measure` (length L1) or its
`measure_2` (length L2), each against an idle of exactly its own length, in one job, every circuit twice (copy 0 / copy 1), shuffled.
Arm S: untouched qubits on the equator, N = 4, 8, read X and Y (the kick). Arm B: measured qubits prepared in 1, N = 4, 8, final
readout of every measured qubit with `measure` (the measured qubit's loss at equal elapsed time). Dry-run by default."""
import os, sys, json, random, time, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit.circuit import Instruction
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
H = os.path.expanduser('~/bones-forge'); HERE = os.path.dirname(os.path.abspath(__file__))
ROLES = H + '/heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_%s.json'
CHIP = sys.argv[1] if len(sys.argv) > 1 else 'ibm_kingston'; FIRE = '--fire' in sys.argv
SHOTS = int(sys.argv[sys.argv.index('--shots') + 1]) if '--shots' in sys.argv else 1000
KEYFILE, INSTANCE = '.ibm_quantum_key', 'open-instance'; SEED = 20261023; NS = (4, 8)
M2 = Instruction('measure_2', 1, 1, [])
def service():
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)
def rot(qc, q, basis):
    if basis == 'X': qc.h(q)
    elif basis == 'Y': qc.sdg(q); qc.h(q)
def circ(n, E, S, mode, N, length, L_dt, basis, arm):
    """mode 'D' idles E for L_dt per step; 'M' reads E with measure (length 'L1') or measure_2 ('L2') per step."""
    qc = QuantumCircuit(n); regs = [ClassicalRegister(len(E), 'm%d' % k) for k in range(N)]
    for r in regs: qc.add_register(r)
    fin = ClassicalRegister(len(S), 'fin'); qc.add_register(fin); efin = None
    if arm == 'B': efin = ClassicalRegister(len(E), 'efin'); qc.add_register(efin); qc.x(E)
    for q in S: qc.sx(q)
    qc.barrier()
    for k in range(N):
        if mode == 'M':
            for j, q in enumerate(E):
                if length == 'L1': qc.measure(q, regs[k][j])
                else: qc.append(M2, [q], [regs[k][j]])
        else:
            for q in E: qc.delay(L_dt, q, unit='dt')
        qc.barrier()
    for i, q in enumerate(S): rot(qc, q, basis); qc.measure(q, fin[i])
    if arm == 'B': qc.measure(E, efin)
    return qc
def build(b):
    n = b.num_qubits; r = json.load(open(ROLES % CHIP)); E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
    L = {'L1': int(round(b.target['measure'][(E[0],)].duration / b.target.dt)), 'L2': int(round(b.target['measure_2'][(E[0],)].duration / b.target.dt))}
    pubs, meta = [], []
    for N in NS:
        for length in ('L1', 'L2'):
            for mode in ('D', 'M'):
                for basis in ('X', 'Y'): pubs.append(circ(n, E, S, mode, N, length, L[length], basis, 'S')); meta.append(dict(arm='S', N=N, length=length, mode=mode, basis=basis))
                pubs.append(circ(n, E, S, mode, N, length, L[length], 'X', 'B')); meta.append(dict(arm='B', N=N, length=length, mode=mode, basis='X'))
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED).shuffle(order)
    return [pubs[i] for i in order], [meta[i] for i in order], E, S, L
if __name__ == '__main__':
    s = service(); b = s.backend(CHIP); pubs, meta, E, S, L = build(b)
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    names = set(op for c in tq for op in c.count_ops()); assert 'measure_2' in names, 'measure_2 did not survive transpilation'
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: L1 {L["L1"]} dt ({L["L1"]*b.target.dt*1e6:.3f} us), L2 {L["L2"]} dt ({L["L2"]*b.target.dt*1e6:.3f} us); {len(tq)} circuits, E {len(E)}, S {len(S)}; '
          f'{dur.min()*1e6:.1f}-{dur.max()*1e6:.1f} us; ESTIMATED {est:.0f} s; remaining {rem}; queue {b.status().pending_jobs}; ops {sorted(names)}')
    if not FIRE: sys.exit(0)
    if rem is not None and rem < est + 3: print('STOP: not enough seconds'); sys.exit(1)
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(f'{HERE}/backend_properties_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = {'chip': CHIP, 'job_id': job.job_id(), 'submitted_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'meta': meta, 'E': E, 'S': S,
           'L_dt': L, 'dt': b.target.dt, 'shots': SHOTS, 'keyfile': KEYFILE, 'instance': INSTANCE, 'plan_sha256': open(HERE + '/PLAN_LENGTH.sha256').read().split()[0]}
    json.dump(rec, open(f'{HERE}/job_{CHIP}.json', 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
