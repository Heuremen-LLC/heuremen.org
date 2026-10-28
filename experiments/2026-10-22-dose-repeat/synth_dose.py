"""synth_dose.py -- plant per-readout kicks in fake S/C-arm data in the control-run raw layout and read them with dose_read.py.
Worlds: none (no kick anywhere), proportional (kick scales with the number read: 39 -> 2 deg, 10 -> 0.51, 1 -> 0.05), global (2 deg for any K)."""
import json, os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import dose_read as DR
def pack(bitsmat):
    n, q = bitsmat.shape; nby = (q + 7) // 8; out = np.zeros((n, nby), np.uint8)
    for k in range(q): out[:, nby - 1 - k // 8] |= (bitsmat[:, k].astype(np.uint8) << (k % 8))
    return out
def world(name, nS=117, shots=1000, seed=5):
    rng = np.random.default_rng(seed); base = rng.uniform(-np.pi, np.pi, nS); amp = 0.8
    kick = {'none': {39: 0, 10: 0, 1: 0}, 'proportional': {39: 2.0, 10: 2.0 * 10 / 39, 1: 2.0 / 39}, 'global': {39: 2.0, 10: 2.0, 1: 2.0}}[name]
    meta = []; raw = {}
    def add(arm, mode, K, st):
        for c in (0, 1):
            i = len(meta); meta.append(dict(arm=arm, mode=mode, N=8, K=K, setting=st, copy=c))
            ph = base + (np.radians(kick[K]) * 8 if mode == 'M' else 0) + (0.02 * rng.normal(size=nS) if c else 0)
            ev = amp * (np.cos(ph) if st == 'XX' else np.sin(ph)); p1 = (1 - ev) / 2
            raw['fin_%d' % i] = pack((rng.random((shots, nS)) < p1).astype(int))
    for mode in ('D', 'M'):
        for st in ('XX', 'YY'): add('S', mode, 39, st)
    for K in (1, 10):
        for mode in ('D', 'M'):
            for st in ('XX', 'YY'): add('C', mode, K, st)
    return dict(chip='synth_' + name, job_id='SYNTH', meta=meta, S=list(range(nS))), raw
if __name__ == '__main__':
    for w in ('none', 'proportional', 'global'):
        rec, raw = world(w); print('### planted world:', w); DR.analyse(rec, raw, sys.stdout)
