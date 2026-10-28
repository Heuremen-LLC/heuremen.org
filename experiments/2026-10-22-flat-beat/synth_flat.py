"""Synthetic test for analyze_flat.py. python synth_flat.py OUTDIR MODE SEED, MODE in track (beat = calibrated zz), flat (beat = 12.5 kHz
whatever the calibration), none (no beat). Writes job/raw/props in the real schema. Physics as synth_zz.py: A dephases with its T2,
B in + splits shots into branches, B relaxes with its T1, B readout flips with its readout error."""
import json, os, sys
import numpy as np
OUT, MODE, SEED = sys.argv[1], sys.argv[2], int(sys.argv[3]); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(SEED)
m = 20; SHOTS = 2000; TAUS = [2 + 4 * k for k in range(41)]; ARMS = ('Z', 'P1', 'P2')
trip = [(3 * j, 3 * j + 1, 3 * j + 2) for j in range(m)]; n = 3 * m
edges = sorted({(3 * j, 3 * j + 1) for j in range(m)} | {(3 * j + 1, 3 * j + 2) for j in range(m)})
zz = {e: rng.uniform(10e3, 57e3) for e in edges}
beat = {e: (zz[e] if MODE == 'track' else 12.5e3 if MODE == 'flat' else 0.0) for e in edges}
T1 = rng.uniform(80, 250, n); T2 = rng.uniform(60, 180, n); ro = rng.uniform(0.005, 0.03, n)
json.dump(dict(qubits=[[dict(name='T1', value=T1[q]), dict(name='T2', value=T2[q]), dict(name='readout_error', value=ro[q])] for q in range(n)],
               general=[dict(name='zz_%d%d' % e, value=zz[e] / 1e9) for e in edges]), open(os.path.join(OUT, 'props_flat_ibm_phoenix.json'), 'w'))
def pack(b):
    nb = (b.shape[1] + 7) // 8; pad = np.zeros((b.shape[0], nb * 8), np.uint8); pad[:, :b.shape[1]] = b
    return np.packbits(pad[:, ::-1], axis=1, bitorder='big')
meta, out, i = [], {}, 0
for c in (0, 1):
    for arm in ARMS:
        for tau in TAUS:
            for st in 'XY':
                t = tau * 1e-6; A = np.zeros((SHOTS, m), np.uint8); r1 = np.zeros_like(A); r2 = np.zeros_like(A)
                for j, (b1, a, b2) in enumerate(trip):
                    phi = rng.normal(0, 1, SHOTS) * np.sqrt(2 * t / (T2[a] * 1e-6))
                    for b, on, reg in ((b1, arm == 'P1', r1), (b2, arm == 'P2', r2)):
                        s1 = (rng.random(SHOTS) < 0.5) if on else np.zeros(SHOTS, bool); td = rng.exponential(T1[b] * 1e-6, SHOTS)
                        phi = phi + 2 * np.pi * beat[(min(a, b), max(a, b))] * np.where(s1, np.minimum(td, t), 0.0)
                        reg[:, j] = (s1 & (td > t)) ^ (rng.random(SHOTS) < ro[b])
                    phi = phi + (0 if st == 'X' else -np.pi / 2); A[:, j] = rng.random(SHOTS) < (1 - np.cos(phi)) / 2
                out['a_%d' % i] = pack(A); out['b1_%d' % i] = pack(r1); out['b2_%d' % i] = pack(r2)
                meta.append(dict(arm=arm, tau_us=tau, tau_dt=0, setting=st, copy=c)); i += 1
np.savez_compressed(os.path.join(OUT, 'raw_flat_ibm_phoenix.npz'), **out)
json.dump(dict(run='flat', chip='ibm_phoenix', n_qubits=n, job_id='SYNTH-%s-%d' % (MODE, SEED), meta=meta, shots=SHOTS, info=dict(triples=trip), edges=edges),
          open(os.path.join(OUT, 'job_flat_ibm_phoenix.json'), 'w'))
print('wrote', OUT, MODE, i)
