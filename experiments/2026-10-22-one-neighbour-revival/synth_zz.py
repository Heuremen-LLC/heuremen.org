"""Synthetic test for analyze_zz.py: writes job/raw/props files in the real schema into a folder and runs nothing on IBM.
python synth_zz.py OUTDIR SCALE SEED   (SCALE = true ZZ scale s; 0 = no ZZ at all, the null). Physics: A dephases with its own T2 envelope and
detuning; B1/B2 in + split each shot into branches; a branch in 1 adds the ZZ frequency s*zz to A until B relaxes (T1); echo inverts the
second half; B readout flips with its readout error."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
OUT, S, SEED = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(SEED)
m = 21; SHOTS = 1000; TAUS = [2 + 5 * k for k in range(41)]; ARMS = ('Z', 'P1', 'P2', 'P12', 'O1', 'EZ', 'EP1')
trip = [(3 * j, 3 * j + 1, 3 * j + 2) for j in range(m)]; n = 3 * m
edges = sorted({(3 * j, 3 * j + 1) for j in range(m)} | {(3 * j + 1, 3 * j + 2) for j in range(m)})
zz = {e: rng.uniform(10e3, 50e3) for e in edges}
T1 = rng.uniform(80, 250, n); T2 = rng.uniform(60, 180, n); ro = rng.uniform(0.005, 0.03, n); det = rng.uniform(-20e3, 20e3, n)
props = dict(qubits=[[dict(name='T1', value=T1[q]), dict(name='T2', value=T2[q]), dict(name='readout_error', value=ro[q])] for q in range(n)],
             general=[dict(name='zz_%d%d' % e, value=zz[e] / 1e9) for e in edges])
json.dump(props, open(os.path.join(OUT, 'props_zz_ibm_phoenix.json'), 'w'))
meta = []; out = {}
def pack(b):   # (shots, m) bits -> qiskit BitArray layout
    nb = (b.shape[1] + 7) // 8; pad = np.zeros((b.shape[0], nb * 8), np.uint8); pad[:, :b.shape[1]] = b
    return np.packbits(pad[:, ::-1], axis=1, bitorder='big')
i = 0
for c in (0, 1):
    for arm in ARMS:
        for tau in TAUS:
            for st in 'XY':
                t = tau * 1e-6; A = np.zeros((SHOTS, m), np.uint8); b1r = np.zeros_like(A); b2r = np.zeros_like(A)
                for j, (b1, a, b2) in enumerate(trip):
                    echo = arm in ('EZ', 'EP1')
                    def branch(on, start1, Tq):
                        st1 = (rng.random(SHOTS) < 0.5) if on == '+' else (np.ones(SHOTS, bool) if on == '1' else np.zeros(SHOTS, bool))
                        td = rng.exponential(Tq * 1e-6, SHOTS); t_in1 = np.where(st1, np.minimum(td, t), 0.0)
                        if echo:   # time in 1 during first half minus during second half
                            h = t / 2; first = np.where(st1, np.minimum(td, h), 0.0); t_in1 = first - (t_in1 - first)
                        end1 = st1 & (td > t); return t_in1, end1
                    s1 = '+' if arm in ('P1', 'P12', 'EP1') else '1' if arm == 'O1' else '0'; s2 = '+' if arm in ('P2', 'P12') else '0'
                    ti1, e1 = branch(s1, 0, T1[b1]); ti2, e2 = branch(s2, 0, T1[b2])
                    k1 = zz[(min(a, b1), max(a, b1))]; k2 = zz[(min(a, b2), max(a, b2))]
                    phi = 2 * np.pi * (S * k1 * ti1 + S * k2 * ti2) + (0 if echo else 2 * np.pi * det[a] * t) + rng.normal(0, 1, SHOTS) * np.sqrt(2 * t / (T2[a] * 1e-6))
                    phi = phi + (0 if st == 'X' else -np.pi / 2)
                    A[:, j] = (rng.random(SHOTS) < (1 - np.cos(phi)) / 2)
                    b1r[:, j] = e1 ^ (rng.random(SHOTS) < ro[b1]); b2r[:, j] = e2 ^ (rng.random(SHOTS) < ro[b2])
                out['a_%d' % i] = pack(A); out['b1_%d' % i] = pack(b1r); out['b2_%d' % i] = pack(b2r)
                meta.append(dict(arm=arm, tau_us=tau, tau_dt=0, setting=st, copy=c)); i += 1
np.savez_compressed(os.path.join(OUT, 'raw_zz_ibm_phoenix.npz'), **out)
json.dump(dict(run='zz', chip='ibm_phoenix', n_qubits=n, job_id='SYNTH-%s-%d' % (S, SEED), meta=meta, shots=SHOTS, info=dict(triples=trip), edges=edges),
          open(os.path.join(OUT, 'job_zz_ibm_phoenix.json'), 'w'))
print('wrote', OUT, 'scale', S, 'pubs', i)
