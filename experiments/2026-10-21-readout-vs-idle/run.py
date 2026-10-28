"""Deciding run (PREREG.md, sha256 in PREREG.sha256). Builds the 30 pubs per chip, dry-runs them, and with --fire submits
one job per chip in order fez, marrakesh, kingston on the Darkwave Disco instance until the remaining seconds run out.
Tokens and CRNs are read from files and never printed."""
import os, sys, json, time, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
H = os.path.expanduser('~/bones-forge'); HERE = os.path.dirname(os.path.abspath(__file__))
ROLES = H + '/heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_%s.json'
CHIPS = sys.argv[sys.argv.index('--chips')+1].split(',') if '--chips' in sys.argv else ['ibm_fez', 'ibm_marrakesh', 'ibm_kingston']; NS = [0, 2, 4, 8, 16]; MODES = ['D', 'M0', 'M1']; QUADS = ['X', 'Y']; SHOTS = 2000
def service():
    tok = open(H + '/.ibm_quantum_key_3').read().strip()
    s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == 'Darkwave Disco'][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)
def build(b, chip):
    r = json.load(open(ROLES % chip)); n = b.num_qubits; E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
    L_dt = int(round(b.target['measure'][(E[0],)].duration / b.target.dt))
    pubs, meta = [], []
    for mode in MODES:
        for N in NS:
            for quad in QUADS:
                mid = ClassicalRegister(len(E), 'mid'); fin = ClassicalRegister(len(S), 'fin')
                qc = QuantumCircuit(n); qc.add_register(mid); qc.add_register(fin)
                if mode == 'M1': qc.x(E)
                qc.sx(S); qc.barrier()
                for _ in range(N):
                    if mode == 'D':
                        for q in E: qc.delay(L_dt, q, unit='dt')
                    else:
                        qc.measure(E, mid)
                    qc.barrier()
                if quad == 'Y': qc.rz(np.pi / 2, S)
                qc.sx(S); qc.measure(S, fin)
                pubs.append(qc); meta.append({'mode': mode, 'N': N, 'quad': quad})
    return pubs, meta, E, S, L_dt
def remaining(s):
    u = s.usage(); return (u.get('usage_remaining_seconds') or 0), u
if __name__ == '__main__':
    fire = '--fire' in sys.argv; s = service(); rem, _ = remaining(s); print('remaining seconds at start:', rem)
    for chip in CHIPS:
        b = s.backend(chip); pubs, meta, E, S, L_dt = build(b, chip)
        tq = transpile(pubs, b, optimization_level=0, initial_layout=list(range(b.num_qubits)), scheduling_method='alap')
        dur = max(c.estimate_duration(b.target, 's') for c in tq)
        print('%s: %d pubs, emitters %d, spectators %d, L=%d dt, longest circuit %.1f us' % (chip, len(tq), len(E), len(S), L_dt, dur * 1e6))
        if not fire: continue
        rem, _ = remaining(s)
        if rem < 10: print('DRY: %.0f s left, stopping before %s' % (rem, chip)); break
        job = SamplerV2(mode=b).run(tq, shots=SHOTS)
        rec = {'chip': chip, 'job_id': job.job_id(), 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'meta': meta, 'E': E, 'S': S, 'L_dt': L_dt, 'shots': SHOTS, 'prereg_sha256': open(HERE + '/PREREG.sha256').read().split()[0]}
        json.dump(rec, open(HERE + '/job_%s.json' % chip, 'w')); print('SUBMITTED', chip, job.job_id(), flush=True)
        res = job.result()
        out = {}
        for i, pr in enumerate(res):
            out['fin_%d' % i] = pr.data.fin.array; out['mid_%d' % i] = pr.data.mid.array if hasattr(pr.data, 'mid') else np.zeros(1)
        np.savez_compressed(HERE + '/raw_%s.npz' % chip, **out)
        m = job.metrics(); used = (m.get('usage') or {}).get('quantum_seconds')
        rem, _ = remaining(s); print('DONE %s: quantum_seconds %s, remaining now %s' % (chip, used, rem), flush=True)
