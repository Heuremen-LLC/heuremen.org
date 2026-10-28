"""Entanglement run (PREREG.md, sha256 in PREREG.sha256). 25 pubs per chip. Dry-run by default; --fire submits.
--acct open4 (free open-plan key_4) or --acct disco (paid Darkwave Disco). Tokens and CRNs are read from files, never printed."""
import os, sys, json, time, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
H = os.path.expanduser('~/bones-forge'); HERE = os.path.dirname(os.path.abspath(__file__))
ROLES = H + '/heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_%s.json'
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
CHIP = arg('--chip', 'ibm_marrakesh'); ACCT = arg('--acct', 'open4'); SHOTS = 1200
CONDS = [('D', 0), ('D', 8), ('D', 16), ('M0', 8), ('M0', 16)]; SETTINGS = ['ZZ', 'XX', 'YY', 'XY', 'YX']
def service():
    if ACCT == 'open4': key, name = '/.ibm_quantum_key_4', 'open-instance'
    else: key, name = '/.ibm_quantum_key_3', 'Darkwave Disco'
    tok = open(H + key).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == name][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)
def pairs_for(b, S):
    Sset = set(S); used = set(); P = []
    for a, c in sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}):
        if a in Sset and c in Sset and a not in used and c not in used: P.append((a, c)); used |= {a, c}
    return P
def rot(qc, q, basis):
    if basis == 'X': qc.h(q)
    elif basis == 'Y': qc.sdg(q); qc.h(q)
def build(b):
    r = json.load(open(ROLES % CHIP)); n = b.num_qubits; E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
    P = pairs_for(b, S); L_dt = int(round(b.target['measure'][(E[0],)].duration / b.target.dt))
    pubs, meta = [], []
    for mode, N in CONDS:
        for st in SETTINGS:
            mid = ClassicalRegister(len(E), 'mid'); fin = ClassicalRegister(2 * len(P), 'fin')
            qc = QuantumCircuit(n); qc.add_register(mid); qc.add_register(fin)
            for a, c in P: qc.h(a); qc.cx(a, c)
            qc.barrier()
            for _ in range(N):
                if mode == 'D':
                    for q in E: qc.delay(L_dt, q, unit='dt')
                else: qc.measure(E, mid)
                qc.barrier()
            for a, c in P: rot(qc, a, st[0]); rot(qc, c, st[1])
            for i, (a, c) in enumerate(P): qc.measure(a, fin[2 * i]); qc.measure(c, fin[2 * i + 1])
            pubs.append(qc); meta.append({'mode': mode, 'N': N, 'setting': st})
    return pubs, meta, E, S, P, L_dt
if __name__ == '__main__':
    fire = '--fire' in sys.argv; s = service(); u = s.usage(); rem = u.get('usage_remaining_seconds') or 0
    print('account', ACCT, 'remaining seconds:', rem)
    b = s.backend(CHIP); pubs, meta, E, S, P, L_dt = build(b)
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = [c.estimate_duration(b.target, 's') for c in tq]
    print('%s: %d pubs, emitters %d, spectators %d, pairs %d, L=%d dt, circuits %.1f-%.1f us, shots %d, total shots %d, queue %s'
          % (CHIP, len(tq), len(E), len(S), len(P), L_dt, min(dur) * 1e6, max(dur) * 1e6, SHOTS, SHOTS * len(tq), b.status().pending_jobs))
    if not fire: sys.exit()
    if rem < 9: print('STOP: %.0f s left' % rem); sys.exit()
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = {'chip': CHIP, 'acct': ACCT, 'job_id': job.job_id(), 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'meta': meta, 'E': E, 'S': S,
           'pairs': P, 'L_dt': L_dt, 'shots': SHOTS, 'prereg_sha256': open(HERE + '/PREREG.sha256').read().split()[0]}
    json.dump(rec, open(HERE + '/job_%s.json' % CHIP, 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res): out['fin_%d' % i] = pr.data.fin.array
    np.savez_compressed(HERE + '/raw_%s.npz' % CHIP, **out)
    print('DONE %s: remaining now %s' % (CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
