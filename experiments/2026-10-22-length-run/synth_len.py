"""synth_len.py -- planted worlds for analyze_len.py, before any job. Chips as kingston (L 2.180 / 1.380 us) and marrakesh (2.684 / 1.884),
5000 shots per circuit. Per untouched qubit: coherence 0.6 (10% dead at 0.02), an idle phase drawn per qubit and N, and a kick per readout
with a qubit-to-qubit spread of 0.8 deg. Worlds: NULL (no kick), PER-READOUT (-2.4 deg at both lengths), PROPORTIONAL (-2.4 deg at L1, scaled
by L2/L1), MIXED (a = -1.0 deg + b = -0.6 deg/us). Measured qubit: P(1) idle 0.80, readout loss 0.03 per 4 readouts."""
import json, os, numpy as np, subprocess, sys
rng = np.random.default_rng(5); SH = 5000; NS = (4, 8); OUT = '/tmp/_lensynth'; os.makedirs(OUT, exist_ok=True)
def pack(b):
    n = b.shape[1]; nby = (n + 7) // 8; a = np.zeros((b.shape[0], nby), np.uint8)
    for k in range(n): a[:, nby - 1 - k // 8] |= (b[:, k].astype(np.uint8) << (k % 8))
    return a
def make(name, L1, L2, kick):
    S = list(range(117)); E = list(range(117, 156)); meta = []
    for N in NS:
        for length in ('L1', 'L2'):
            for mode in ('D', 'M'):
                for basis in ('X', 'Y'): meta.append(dict(arm='S', N=N, length=length, mode=mode, basis=basis))
                meta.append(dict(arm='B', N=N, length=length, mode=mode, basis='X'))
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]
    coh = np.where(rng.random(117) < 0.1, 0.02, 0.6); base = {N: rng.uniform(-np.pi, np.pi, 117) for N in NS}; spread = rng.normal(0, 0.8, 117)
    raw = {}
    for i, m in enumerate(meta):
        N, L = m['N'], (L1 if m['length'] == 'L1' else L2)
        if m['arm'] == 'S':
            ph = base[N] + (np.radians(N * (kick(L) + spread)) if m['mode'] == 'M' else 0)
            ex = coh * (np.cos(ph) if m['basis'] == 'X' else np.sin(ph)); p1 = (1 - ex) / 2
            raw[f'fin_{i}'] = pack(rng.random((SH, 117)) < p1)
        else:
            p1 = 0.80 - (0.03 * N / 4 if m['mode'] == 'M' else 0); raw[f'fin_{i}'] = pack(rng.random((SH, 117)) < 0.5); raw[f'efin_{i}'] = pack(rng.random((SH, 39)) < p1)
    json.dump({'chip': name, 'job_id': 'synthetic', 'meta': meta, 'E': E, 'S': S, 'L_dt': {'L1': round(L1 / 0.004), 'L2': round(L2 / 0.004)}, 'dt': 4e-9, 'shots': SH}, open(f'{OUT}/job_{name}.json', 'w'))
    np.savez_compressed(f'{OUT}/raw_{name}.npz', **raw)
worlds = {'NULL': lambda L: 0.0, 'PERREADOUT': lambda L: -2.4, 'PROP': None, 'MIXED': lambda L: -1.0 - 0.6 * L}
names = []
for chip, L1, L2 in (('kingston', 2.180, 1.380), ('marrakesh', 2.684, 1.884)):
    for w, f in worlds.items():
        if w == 'PROP': f = (lambda L1_: (lambda L: -2.4 * L / L1_))(L1)
        nm = f'synth_{w}_{chip}'; make(nm, L1, L2, f); names.append(nm)
print(subprocess.run([sys.executable, 'analyze_len.py', *names, '--dir', OUT], capture_output=True, text=True).stdout)
