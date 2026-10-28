"""Synthetic test of analyze_dd.py. Builds dd.py's real circuit list on FakeFez (offline, nothing submitted), draws bits for every pub in the
schema a fired job writes (job_dd_<chip>.json, raw_dd_<chip>.npz with fin_<i>, props_dd_<chip>.json), runs the analysis.
  python synth_dd.py            -> PLANTED (folder synth_planted) then NULL seed 5 (synth_null_5); both outputs into synth_output.txt
  python synth_dd.py --nulls    -> NULL seeds 11-16 (nulls/), tally of resolved readings into null_floor.txt
  MODE=planted|null SEED=k python synth_dd.py --one   -> one case to stdout (used by the two above)
World (both cases): each spectator has its own coherence c0 in [0.55, 0.95] (6 spectators dead at 0.01), decay 0.4 % per step, its own
starting phase, and its own drift per step -- DIFFERENT drift with the echo on (the echo refocuses detuning), the same in D and M.
PLANTED: kick -2.0 deg per readout with DD off, 0 with DD on (removed); coherence cost 0.010 at every N in BOTH DD settings (not removed);
echo has no effect of its own on the idle arm. NULL: no kick, no cost. Only the 'fin' register is drawn (the mid-circuit m<k> registers are not read)."""
import os, sys, json, subprocess, warnings; warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable
def one(mode, seed):
    sys.argv = sys.argv[:1]; sys.path.insert(0, HERE)
    import dd, analyze_dd
    from qiskit_ibm_runtime.fake_provider import FakeFez
    A = analyze_dd.A; A.RNG = np.random.default_rng(seed + 1000)
    b = FakeFez(); n = b.num_qubits; rng = np.random.default_rng(seed); chip = 'ibm_fez'
    OUT = os.path.join(HERE, 'synth_planted' if mode == 'planted' else ('synth_null_%d' % seed if seed == 5 else os.path.join('nulls', 'synth_null_%d' % seed)))
    os.makedirs(OUT, exist_ok=True)
    pubs, meta, shots, info = dd.build(b, chip); S = info['S']; nS = len(S)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]
    c0 = rng.uniform(0.55, 0.95, nS); c0[rng.choice(nS, 6, replace=False)] = 0.01
    phi0 = rng.uniform(-np.pi, np.pi, nS); drift = {'off': np.radians(rng.normal(0, 8, nS)), 'on': np.radians(rng.normal(0, 1, nS))}
    KICK = {'off': np.radians(-2.0), 'on': 0.0} if mode == 'planted' else {'off': 0.0, 'on': 0.0}; COST = 0.010 if mode == 'planted' else 0.0
    z = {}
    for i, m in enumerate(meta):
        N = m['N']; c = c0 * np.exp(-0.004 * N) - (COST if m['mode'] == 'M' else 0.0); phi = phi0 + drift[m['dd']] * N + (KICK[m['dd']] * N if m['mode'] == 'M' else 0.0)
        e = c * (np.cos(phi) if m['setting'] == 'XX' else np.sin(phi)); p = np.zeros(n); p[S] = (1 - e) / 2
        z['fin_%d' % i] = (rng.random((shots, n)) < p).astype(np.int8)
    props = b.properties().to_dict(); json.dump(props, open(os.path.join(OUT, 'props_dd_%s.json' % chip), 'w'), default=str)
    rec = dict(run='dd', chip=chip, n_qubits=n, job_id='SYNTH_%s_%d' % (mode.upper(), seed), meta=meta, order=list(range(len(meta))), shots=shots, info=info, plan_sha256='SYNTH')
    json.dump(rec, open(os.path.join(OUT, 'job_dd_%s.json' % chip), 'w'), default=str); np.savez_compressed(os.path.join(OUT, 'raw_dd_%s.npz' % chip), **z)
    if mode == 'planted':
        print('PLANTED (seed %d): kick -2.0 deg/readout DD off, 0 DD on -> expect R1 BELOW ~-2, R2 CONSISTENT, R3 ABOVE ~+2, R3s BELOW ~-2;'
              ' coherence cost 0.010 both settings (minus a small hypot bias) -> R4, R5 ABOVE ~0.01, R6 CONSISTENT; R7 CONSISTENT; R0 CONSISTENT; R8: 6 dead' % seed)
    else: print('NULL (seed %d): nothing planted; every reading has truth 0' % seed)
    analyze_dd.run_all(OUT)
if __name__ == '__main__':
    if '--one' in sys.argv: one(os.environ.get('MODE', 'null'), int(os.environ.get('SEED', '5'))); sys.exit()
    def sub(mode, seed):
        return subprocess.run([PY, os.path.abspath(__file__), '--one'], env=dict(os.environ, MODE=mode, SEED=str(seed)), capture_output=True, text=True, check=True).stdout
    if '--nulls' in sys.argv:
        os.makedirs(os.path.join(HERE, 'nulls'), exist_ok=True); tot = np.zeros(4, int); lines = []
        for sd in range(11, 17):
            out = sub('null', sd); open(os.path.join(HERE, 'nulls', 'null_%d.txt' % sd), 'w').write(out)
            t = [l for l in out.splitlines() if l.startswith('TALLY')][0]; lines.append('seed %d: %s' % (sd, t))
            k = [int(x.split()[-1]) for x in t.split(':', 1)[1].split(',')]; tot += np.array(k)
        s = ('NULL FLOOR, 6 seeds (11-16), every reading truth 0, all/wired/unwired subsets counted separately\n' + '\n'.join(lines) +
             '\nTOTAL readings %d: ABOVE %d, BELOW %d -> resolved %d (%.1f %%); CONSISTENT %d; CANNOT TELL %d\n' % (tot.sum(), tot[0], tot[1], tot[0] + tot[1], 100 * (tot[0] + tot[1]) / tot.sum(), tot[2], tot[3]))
        open(os.path.join(HERE, 'null_floor.txt'), 'w').write(s); print(s)
    else:
        out = sub('planted', 5) + '\n' + '#' * 100 + '\n' + sub('null', 5); open(os.path.join(HERE, 'synth_output.txt'), 'w').write(out); print(out)
