"""run_threebox.py -- three-box weak-value run (PLAN_THREEBOX.md). Usage: run_threebox.py CHIP [--fire] [--shots N] [--key KEYFILE] [--instance NAME]
Dry-run by default: builds, transpiles, prints layout and estimated seconds, submits nothing. --fire submits one job (Darkwave Disco),
saves backend properties at firing, then computes IBM's noise-model prediction for every transpiled circuit from those properties BEFORE
the job's result is read, and saves it (model_<chip>.json)."""
import os, sys, json, random, warnings, time; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from tb_common import EPS, BOXES, VARIANTS, circuit, keys
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
CHIP = sys.argv[1]; FIRE = '--fire' in sys.argv; SHOTS = int(sys.argv[sys.argv.index('--shots') + 1]) if '--shots' in sys.argv else 4000
SEED = 20261022; OVERHEAD = 290e-6; FEZ_LAYOUT = [142, 143, 144]          # s0, s1, probe: the qubits of the 11 Oct fez jobs
KEYFILE = sys.argv[sys.argv.index('--key') + 1] if '--key' in sys.argv else '.ibm_quantum_key_3'
INSTANCE = sys.argv[sys.argv.index('--instance') + 1] if '--instance' in sys.argv else 'Darkwave Disco'
def service():
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)
def choose_layout(b):
    """fixed rule: fez uses the 11 Oct qubits; elsewhere the chain a-m-c (edges a-m, m-c) with the lowest sum of readout errors of a, m, c
    plus the two two-qubit gate errors, from the properties at firing; s1 = middle m, s0 = lower-index end, probe = other end."""
    if b.name == 'ibm_fez': return FEZ_LAYOUT
    t = b.target; edges = {tuple(sorted(e)) for e in b.coupling_map.get_edges()}; nb = {}
    for a, c in edges: nb.setdefault(a, set()).add(c); nb.setdefault(c, set()).add(a)
    g2 = [n for n in ('cz', 'ecr', 'cx') if n in t.operation_names][0]
    def ro(q):
        p = t['measure'][(q,)]; return p.error if p and p.error is not None else 1.0
    def ge(a, c):
        for x in ((a, c), (c, a)):
            if x in t[g2] and t[g2][x] and t[g2][x].error is not None: return t[g2][x].error
        return 1.0
    best = None
    for m, ns in nb.items():
        ns = sorted(ns)
        for i in range(len(ns)):
            for j in range(i + 1, len(ns)):
                a, c = ns[i], ns[j]; cost = ro(a) + ro(m) + ro(c) + ge(a, m) + ge(m, c)
                if best is None or cost < best[0]: best = (cost, [a, m, c])
    a, m, c = best[1]; return [min(a, c), m, max(a, c)]
def build(b):
    lay = choose_layout(b); K = keys(); items = [(kk, cp) for kk in K for cp in (0, 1)]; random.Random(SEED).shuffle(items)
    logical = {kk: circuit(*kk) for kk in K}
    tq = transpile([logical[kk] for kk, _ in items], backend=b, initial_layout=lay, optimization_level=1, seed_transpiler=SEED)
    return lay, items, logical, tq
def ideal_xp(qc):
    p = Statevector(qc.remove_final_measurements(inplace=False)).probabilities([0, 1, 2]); a, c = p[0b000], p[0b100]; return float((a - c) / (a + c))
def model_xp(b, tq, items):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel
    sim = AerSimulator(method='density_matrix', noise_model=NoiseModel.from_backend(b))
    out = {}
    for (kk, cp), c in zip(items, tq):
        if cp: continue
        c2 = c.copy(); c2.remove_final_measurements(); act = sorted({c2.find_bit(q).index for ins in c2.data for q in ins.qubits})
        meas = [c.find_bit(ins.qubits[0]).index for ins in c.data if ins.operation.name == 'measure']   # physical s0, s1, probe in clbit order
        c4 = c2.copy(); c4.save_density_matrix(qubits=meas); dm = sim.run(c4).result().data()['density_matrix']
        pr = np.real(np.diag(np.asarray(dm))); a, z = pr[0b000], pr[0b100]; out[f'{kk[0]}|{kk[1]}|{kk[2]:.2f}'] = float((a - z) / (a + z))
    return out
if __name__ == '__main__':
    s = service(); b = s.backend(CHIP); lay, items, logical, tq = build(b)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + OVERHEAD))
    print(f'{CHIP}: layout s0,s1,probe = {lay}; {len(tq)} circuits ({len(keys())} x 2 copies); shots/circuit {SHOTS}; circuit {dur.min()*1e6:.1f}-{dur.max()*1e6:.1f} us; '
          f'ESTIMATED {est:.0f} s; 2q gates per circuit {min(sum(v for k_, v in c.count_ops().items() if k_ in ("cz","ecr","cx")) for c in tq)}-'
          f'{max(sum(v for k_, v in c.count_ops().items() if k_ in ("cz","ecr","cx")) for c in tq)}; queue {b.status().pending_jobs}; remaining s {s.usage().get("usage_remaining_seconds")}')
    if not FIRE: sys.exit(0)
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(f'{HERE}/backend_properties_{CHIP}.json', 'w'), default=str)
    smp = SamplerV2(mode=b); job = smp.run(tq, shots=SHOTS); jid = job.job_id()
    meta = {'job': jid, 'chip': CHIP, 'keyfile': KEYFILE, 'instance': INSTANCE, 'layout': lay, 'shots': SHOTS, 'items': [[list(kk), cp] for kk, cp in items], 'fired_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'plan_sha256': open(HERE + '/PLAN_THREEBOX.sha256').read().split()[0]}
    json.dump(meta, open(f'{HERE}/job_{CHIP}.json', 'w'), indent=1); print('submitted', jid, flush=True)
    model = {'chip': CHIP, 'model': model_xp(b, tq, items), 'ideal': {f'{k[0]}|{k[1]}|{k[2]:.2f}': ideal_xp(c) for k, c in logical.items()},
             'computed_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'job_status_when_computed': str(job.status())}
    json.dump(model, open(f'{HERE}/model_{CHIP}.json', 'w'), indent=1); print('IBM noise-model prediction saved before reading results; job status', model['job_status_when_computed'])
