"""Synthetic test of analyze_reach.py. Builds the real reach circuits list on FakeFez, FakeKingston, FakeMarrakesh (offline; same layout rule),
plants a known kick, writes the files a fired job writes (job_reach_<chip>.json, raw_reach_<chip>.npz bit-packed as SamplerV2 returns them,
props_reach_<chip>.json), and runs the analysis. Scenarios: null (no kick), step (kick to d=5, 0 beyond), smooth (decays with d), and 6 null seeds.
Realism planted: per-qubit coherence 0.45-0.80 (readout arm 3 % lower), 3 dead clocks per chip (coherence 0.02), a static phase per qubit
(-60..60 deg, both arms), a per-unit spread of the kick (x exp(0.5 g), g normal: the real d=3-5 scatter is large). Copies identical in truth.
python synth_reach.py  -> writes synth_reach/<scenario>/ and synth_reach/SYNTH_RESULTS.txt"""
import os, sys, json, io, contextlib, warnings; warnings.filterwarnings('ignore'); sys.argv = sys.argv[:1]
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import reach, analyze_reach as AR
from qiskit_ibm_runtime.fake_provider import FakeFez, FakeKingston, FakeMarrakesh
CHIPS = (('ibm_fez', FakeFez), ('ibm_kingston', FakeKingston), ('ibm_marrakesh', FakeMarrakesh))
BUILT = {}
def planted(kind, d):
    d = np.asarray(d, float)
    if kind == 'null': return 0 * d
    if kind == 'step': return np.where(d <= 5, -2.0, 0.0)
    if kind == 'smooth': return -3.0 * np.exp(-(d - 1) / 2.5)
def pack(b):
    nb = (b.shape[1] + 7) // 8; full = np.zeros((b.shape[0], nb * 8), np.uint8); full[:, :b.shape[1]] = b
    return np.packbits(full[:, ::-1], axis=1, bitorder='big')
def make(kind, seed, out):
    rng = np.random.default_rng(seed); os.makedirs(out, exist_ok=True); truth = {}
    for chip, F in CHIPS:
        if chip not in BUILT:
            b = F(); p, m, s, info = reach.build(b, chip, verbose=False); BUILT[chip] = (b.num_qubits, m, s, info, json.loads(json.dumps(b.properties().to_dict(), default=str)))
        n, meta0, shots, info = BUILT[chip][:4]; props = BUILT[chip][4]
        meta = [dict(m, copy=c) for c in (0, 1) for m in meta0]
        co = rng.uniform(0.45, 0.80, n); dead = rng.choice(n, 3, replace=False); co[dead] = 0.02; ph0 = np.radians(rng.uniform(-60, 60, n))
        spread = [np.exp(0.5 * rng.standard_normal(len(l['S']))) for l in info['layouts']]
        kk = [planted(kind, l['d']) * f for l, f in zip(info['layouts'], spread)]   # deg per readout, per unit
        units_d = np.concatenate([l['d'] for l in info['layouts']]); units_alive = ~np.isin(np.concatenate([l['S'] for l in info['layouts']]), dead)
        truth[chip] = dict(size_per_unit=-np.concatenate(kk), d=units_d, alive=units_alive)
        z = {}
        for i, m in enumerate(meta):
            if m['layout'] == -1: S = np.arange(n); c = co; phi = ph0
            else:
                S = np.array(info['layouts'][m['layout']]['S']); c = co[S] * (0.97 if m['mode'] == 'M' else 1.0)
                phi = ph0[S] + (np.radians(kk[m['layout']]) * m['N'] if m['mode'] == 'M' else 0)
            e = c * (np.cos(phi) if m['setting'] == 'XX' else np.sin(phi)); bits = (rng.random((shots, len(S))) < (1 - e) / 2).astype(np.uint8)
            z['fin_%d' % i] = pack(bits)
        json.dump(dict(run='reach', chip=chip, n_qubits=n, job_id='SYNTH-%s-%d' % (kind, seed), meta=meta, info=info, shots=shots, plan_sha256='SYNTH'),
                  open(os.path.join(out, 'job_reach_%s.json' % chip), 'w'), default=str)
        json.dump(props, open(os.path.join(out, 'props_reach_%s.json' % chip), 'w')); np.savez_compressed(os.path.join(out, 'raw_reach_%s.npz' % chip), **z)
    return truth
def truth_bins(t):
    b = np.minimum(t['d'], 8); return {k: float(np.median(t['size_per_unit'][(b == k) & t['alive']])) for k in AR.BINS}
def run(kind, seed, folder):
    truth = make(kind, seed, folder); AR.A.RNG = np.random.default_rng(seed + 1000); buf = io.StringIO()
    with contextlib.redirect_stdout(buf): per, comb = AR.run_all(folder)
    open(os.path.join(folder, 'analysis_output.txt'), 'w').write(buf.getvalue()); return truth, buf.getvalue(), comb
if __name__ == '__main__':
    R = []; P = lambda *a: (print(*a), R.append(' '.join(str(x) for x in a)))
    TOP = os.path.join(HERE, 'synth_reach')
    for kind in ('null', 'step', 'smooth'):
        truth, txt, comb = run(kind, 5, os.path.join(TOP, kind))
        P('=' * 100); P('SCENARIO %s (seed 5). Planted size per bin = median over alive units of the planted per-unit size (deg per readout).' % kind)
        P('  planted kick: ' + {'null': '0 everywhere', 'step': '-2.0 deg for d<=5, 0 beyond (x per-unit spread)', 'smooth': '-3.0 x exp(-(d-1)/2.5) deg (x per-unit spread)'}[kind])
        for line in txt.splitlines():
            if line.startswith('reach ') or 'RR' in line or 'guard' in line or 'combined' in line: P(line)
        for chip in truth:
            tb = truth_bins(truth[chip]); P('  TRUTH %s: ' % chip.replace('ibm_', '') + ' '.join('%s %+.3f' % (AR.blab(k), v) for k, v in tb.items()))
        tb = {k: np.mean([truth_bins(truth[c])[k] for c in truth]) for k in AR.BINS}; cov = []
        for k in AR.BINS:
            c = comb.get('RR1 kick size at %s' % AR.blab(k))
            if c: cov.append('%s truth %+.3f read %+.3f [%+.3f, %+.3f] %s %s' % (AR.blab(k), tb[k], c[0], c[1], c[2], c[3], 'covers' if c[1] <= tb[k] <= c[2] else 'MISSES'))
        P('  COMBINED vs mean planted: ' + ' | '.join(cov))
    P('=' * 100); P('NULL FLOOR: six null seeds, no kick planted, every reading (per chip RR0, RR1 x6, RR2 x5; and the same 12 combined)')
    tot = dict(chip=[0, 0, 0, 0], comb=[0, 0, 0, 0]); EDGES = [0, 0]   # resolved, above, below, total; RR3 edges claimed (per chip, combined)
    for seed in range(11, 17):
        truth, txt, comb = run('null', seed, os.path.join(TOP, 'null_seed_%d' % seed))
        chiptxt = txt.split('ACROSS CHIPS')[0]; lines = [l for l in chiptxt.splitlines() if l.strip().startswith(('RR0', 'RR1', 'RR2'))]
        e_chip = sum(1 for l in chiptxt.splitlines() if 'reach edge' in l and 'edge (unbroken ABOVE from d=3) none' not in l)
        e_comb = int('combined: reach edge (unbroken ABOVE from d=3) none' not in txt); EDGES[0] += e_chip; EDGES[1] += e_comb
        a = sum('-> ABOVE' in l for l in lines); bl = sum('-> BELOW' in l for l in lines); ct = sum('CANNOT TELL' in l for l in lines)
        ca = sum(1 for v in comb.values() if v and v[3] == 'ABOVE'); cb = sum(1 for v in comb.values() if v and v[3] == 'BELOW')
        rr3 = [l.strip() for l in txt.splitlines() if 'RR3' in l]
        P('  seed %d: per-chip readings %d, ABOVE %d, BELOW %d, CANNOT TELL (copies) %d | combined %d, ABOVE %d, BELOW %d | %s'
          % (seed, len(lines), a, bl, ct, len(comb), ca, cb, rr3[-1].split('   (')[0]))
        tot['chip'] = [tot['chip'][0] + a + bl, tot['chip'][1] + a, tot['chip'][2] + bl, tot['chip'][3] + len(lines)]
        tot['comb'] = [tot['comb'][0] + ca + cb, tot['comb'][1] + ca, tot['comb'][2] + cb, tot['comb'][3] + len(comb)]
    for k, v in tot.items(): P('  FLOOR %s: %d of %d readings resolved (%.1f %%; %d above, %d below)' % (k, v[0], v[3], 100 * v[0] / v[3], v[1], v[2]))
    P('  RR3 under no kick: a reach edge was named on %d of 18 chip-runs and %d of 6 combined (the edge needs d=3 ABOVE first)' % tuple(EDGES))
    open(os.path.join(TOP, 'SYNTH_RESULTS.txt'), 'w').write('\n'.join(R) + '\n')
