"""synth_direct.py -- planted worlds for analyze_direct.py on the phoenix triples and their stated zz (zz_flat props): TRACK (shift = stated),
HALF (0.5 x stated, a convention factor), FLAT (13 kHz on every edge). Per shot: A detuning (per qubit, same in all arms), coherence decay
T2 = 60 us, neighbour in 1 decays with T1 = 120 us (the shift stops when it decays; heralding on its final readout), 1000 shots, two copies."""
import json, os, sys, numpy as np, subprocess
rng = np.random.default_rng(3); OUT = '/tmp/_dsynth'; os.makedirs(OUT, exist_ok=True)
A1 = os.path.dirname(os.path.abspath(__file__)) + '/..'
src = json.load(open(A1 + '/zz_flat/job_flat_ibm_phoenix.json')); props = json.load(open(A1 + '/zz_flat/props_flat_ibm_phoenix.json'))
sys.argv = ['x']; import analyze_direct as AD
trip = [tuple(t) for t in src['info']['triples']]; edges = {tuple(e) for e in src['edges']}; zz = AD.zz_table(props, edges); m = len(trip)
TAUS = [1 + 2 * k for k in range(20)]; SH = 1000; T2, T1 = 60e-6, 120e-6
def pack(bits):
    n = bits.shape[1]; padded = np.zeros((bits.shape[0], ((n + 7) // 8) * 8), np.uint8); padded[:, :n] = bits
    return np.packbits(padded[:, ::-1][:, -((n + 7) // 8) * 8:], axis=1, bitorder='big') if False else np.packbits(np.concatenate([np.zeros((bits.shape[0], (8 - n % 8) % 8), np.uint8), bits[:, ::-1]], 1), axis=1, bitorder='big')
def world(name, shift):
    det = rng.uniform(-60e3, 60e3, m); meta = []; raw = {}
    for arm in ('Z', 'O1', 'O2'):
        for t in TAUS:
            for st in 'XY': meta.append(dict(arm=arm, tau_us=t, setting=st))
    meta = [dict(x, copy=c) for c in (0, 1) for x in meta]
    for i, mt in enumerate(meta):
        t = mt['tau_us'] * 1e-6; phi = 2 * np.pi * det * t; b1 = np.zeros((SH, m), np.uint8); b2 = np.zeros((SH, m), np.uint8)
        ph = np.tile(phi, (SH, 1))
        for arm, breg, col in (('O1', b1, 0), ('O2', b2, 2)):
            if mt['arm'] != arm: continue
            f = np.array([shift(zz.get((min(tr[1], tr[col]), max(tr[1], tr[col])), 0.0)) for tr in trip])
            tj = rng.exponential(T1, (SH, m)); ph = ph + 2 * np.pi * f * np.minimum(tj, t); breg[:] = (tj > t)
        c = np.exp(-t / T2); ex = c * (np.cos(ph) if mt['setting'] == 'X' else np.sin(ph))
        a = (rng.random((SH, m)) < (1 - ex) / 2).astype(np.uint8)
        raw[f'a_{i}'] = pack(a); raw[f'b1_{i}'] = pack(b1); raw[f'b2_{i}'] = pack(b2)
    json.dump(dict(chip=name, job_id='synthetic', meta=meta, shots=SH, triples=trip, edges=sorted(edges)), open(f'{OUT}/job_direct_{name}.json', 'w'))
    np.savez_compressed(f'{OUT}/raw_direct_{name}.npz', **raw); json.dump(props, open(f'{OUT}/props_direct_{name}.json', 'w'))
for nm, f in (('synthTRACK', lambda s: s), ('synthHALF', lambda s: 0.5 * s), ('synthFLAT', lambda s: 13e3)): world(nm, f)
print(subprocess.run([sys.executable, 'analyze_direct.py', 'synthTRACK', 'synthHALF', 'synthFLAT', '--dir', OUT], capture_output=True, text=True).stdout)
