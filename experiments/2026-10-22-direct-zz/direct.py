"""direct.py -- direct per-edge zz (PLAN_DIRECT.md). Usage: direct.py CHIP [--fire] [--shots N]. Dry run by default.
Triples (B1, A, B2) fixed in advance: ibm_phoenix = the 20 triples of the flat-beat job (zz_flat/job_flat_ibm_phoenix.json); ibm_fez = the
33 triples of the one-neighbour job (zz_revival/job_zz_ibm_fez.json). A: SX, wait tau, read X or Y. Arms: Z (B1, B2 in 0), O1 (B1 in 1),
O2 (B2 in 1). tau = 1-39 us in 2 us steps. Every circuit twice (copy 0 / copy 1), shuffled together. Open instance of .ibm_quantum_key."""
import os, sys, json, random, time, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge'); A1 = os.path.join(HERE, '..')
SRC = {'ibm_phoenix': A1 + '/zz_flat/job_flat_ibm_phoenix.json', 'ibm_fez': A1 + '/zz_revival/job_zz_ibm_fez.json'}
ARMS = ('Z', 'O1', 'O2'); TAU_US = [1 + 2 * k for k in range(20)]; SEED = 20261025
CHIP = sys.argv[1]; FIRE = '--fire' in sys.argv; SHOTS = int(sys.argv[sys.argv.index('--shots') + 1]) if '--shots' in sys.argv else 1000
KEYFILE, INSTANCE = '.ibm_quantum_key', 'open-instance'
def service():
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)
def tau_dt(b, us): g = 16; return int(round(us * 1e-6 / b.target.dt / g)) * g
def rot(qc, q, st):
    if st == 'X': qc.h(q)
    else: qc.sdg(q); qc.h(q)
def build(b):
    trip = [tuple(t) for t in json.load(open(SRC[CHIP]))['info']['triples']]; n = b.num_qubits
    A = [t[1] for t in trip]; B1 = [t[0] for t in trip]; B2 = [t[2] for t in trip]; pubs, meta = [], []
    for arm in ARMS:
        for tus in TAU_US:
            T = tau_dt(b, tus)
            for st in 'XY':
                qc = QuantumCircuit(n); ra, r1, r2 = (ClassicalRegister(len(A), x) for x in ('a', 'b1', 'b2'))
                for r in (ra, r1, r2): qc.add_register(r)
                if arm == 'O1': qc.x(B1)
                if arm == 'O2': qc.x(B2)
                qc.sx(A); qc.barrier(); qc.delay(T, A + B1 + B2, unit='dt'); qc.barrier()
                for q in A: rot(qc, q, st)
                qc.measure(A, ra); qc.measure(B1, r1); qc.measure(B2, r2)
                pubs.append(qc); meta.append(dict(arm=arm, tau_us=tus, tau_dt=T, setting=st))
    return pubs, meta, trip
if __name__ == '__main__':
    s = service(); b = s.backend(CHIP); pubs, meta, trip = build(b)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED).shuffle(order); pubs = [pubs[i] for i in order]; meta = [meta[i] for i in order]
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: {len(trip)} triples, {len(tq)} circuits, shots {SHOTS}, ESTIMATED {est:.0f} s, remaining {rem}, queue {b.status().pending_jobs}')
    if not FIRE: sys.exit(0)
    if rem is not None and rem < est + 3: print('STOP: not enough seconds'); sys.exit(1)
    props = b.properties().to_dict(); json.dump(props, open(f'{HERE}/props_direct_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), meta=meta, shots=SHOTS, triples=trip,
               edges=sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}), keyfile=KEYFILE, instance=INSTANCE, plan_sha256=open(HERE + '/PLAN_DIRECT.sha256').read().split()[0])
    json.dump(rec, open(f'{HERE}/job_direct_{CHIP}.json', 'w'), default=str); print('SUBMITTED', CHIP, job.job_id(), flush=True)
