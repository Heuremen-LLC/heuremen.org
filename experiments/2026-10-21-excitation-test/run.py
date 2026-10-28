"""Excitation test (PREREG.md, sha256 in PREREG.sha256). Darkwave Disco instance. Dry-run by default; --fire submits."""
import os, sys, json, time, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
H = os.path.expanduser('~/bones-forge'); HERE = os.path.dirname(os.path.abspath(__file__))
ROLES = H + '/heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_%s.json'
CHIP = sys.argv[sys.argv.index('--chip') + 1] if '--chip' in sys.argv else 'ibm_kingston'
MODES = ['D0', 'M0', 'D1', 'M1']; NS = [0, 4, 8]; QUADS = ['X', 'Y']; SHOTS = 1500
def service():
    tok = open(H + '/.ibm_quantum_key_3').read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == 'Darkwave Disco'][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)
def build(b):
    r = json.load(open(ROLES % CHIP)); n = b.num_qubits; E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
    L_dt = int(round(b.target['measure'][(E[0],)].duration / b.target.dt)); pubs, meta = [], []
    for mode in MODES:
        for N in NS:
            for quad in QUADS:
                fin = ClassicalRegister(len(S), 'fin'); regs = [ClassicalRegister(len(E), 'm%d' % k) for k in range(max(N, 1))]
                qc = QuantumCircuit(n); [qc.add_register(x) for x in regs]; qc.add_register(fin)
                if mode in ('D1', 'M1'): qc.x(E)
                qc.sx(S); qc.barrier()
                for k in range(N):
                    if mode.startswith('D'):
                        for q in E: qc.delay(L_dt, q, unit='dt')
                    else: qc.measure(E, regs[k])
                    qc.barrier()
                if quad == 'Y': qc.rz(np.pi / 2, S)
                qc.sx(S); qc.measure(S, fin)
                pubs.append(qc); meta.append({'mode': mode, 'N': N, 'quad': quad})
    return pubs, meta, E, S, L_dt
if __name__ == '__main__':
    fire = '--fire' in sys.argv; s = service(); rem = s.usage().get('usage_remaining_seconds') or 0; print('remaining seconds:', rem)
    b = s.backend(CHIP); pubs, meta, E, S, L_dt = build(b)
    tq = transpile(pubs, b, optimization_level=0, initial_layout=list(range(b.num_qubits)), scheduling_method='alap')
    dur = [c.estimate_duration(b.target, 's') for c in tq]
    print('%s: %d pubs, E %d, S %d, L=%d dt, circuits %.1f-%.1f us, shots %d (total %d), queue %s' % (CHIP, len(tq), len(E), len(S), L_dt, min(dur)*1e6, max(dur)*1e6, SHOTS, SHOTS*len(tq), b.status().pending_jobs))
    if not fire: sys.exit()
    if rem < 12: print('STOP: %.0f s left' % rem); sys.exit()
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = {'chip': CHIP, 'job_id': job.job_id(), 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'meta': meta, 'E': E, 'S': S, 'L_dt': L_dt, 'shots': SHOTS, 'prereg_sha256': open(HERE + '/PREREG.sha256').read().split()[0]}
    json.dump(rec, open(HERE + '/job_%s.json' % CHIP, 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        d = pr.data; out['fin_%d' % i] = d.fin.array
        for k in range(8):
            if hasattr(d, 'm%d' % k) and meta[i]['mode'].startswith('M') and k < meta[i]['N']: out['m%d_%d' % (k, i)] = getattr(d, 'm%d' % k).array
    np.savez_compressed(HERE + '/raw_%s.npz' % CHIP, **out)
    print('DONE %s: remaining now %s' % (CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
