"""ladder.py -- the coherence ladder: T1, T2* (Ramsey), T2 echo (Hahn), CPMG-4, CPMG-16 on the SAME qubits in ONE job (PLAN_LADDER.md).
  ladder.py CHIP                 dry run: build, transpile, verify pulse counts, estimate seconds
  ladder.py CHIP --fire          submit (instance INSTANCE below)
  ladder.py --fetch CHIP         save raw_ladder_CHIP.npz
  ladder.py --synth WORLD        the analysis on planted data (worlds in synth())
  ladder.py --analyze CHIP       the pre-registered readings
Circuits, per measured qubit (one sparse set, no two measured qubits coupled, every other qubit idle in |0>):
  T1      X, wait t, measure                                      P1 = B + A exp(-t/T1)
  RAMSEY  SX, wait t, [RZ(pi/2)], SXdg, measure  (two quadratures) coherence c = |(1-2Pa, 1-2Pb)| = A exp(-t/T2*), detuning-free
  HAHN    SX, wait t/2, X, wait t/2, SXdg (= hahn.py)              P1 = 1/2 + A exp(-t/T2)
  CPMG-N  SX, wait t/2N, [Y, wait t/N] x N (last wait t/2N), SXdg  1 - P1 = 1/2 + A exp(-t/T2(N)); N = 4, 16; barriers keep every pulse
Every circuit twice (copy 0 / copy 1), shuffled together."""
import os, sys, json, time, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
sys.path.insert(0, os.path.join(HERE, '..', 'hahn_t2'))
from hahn import sparse_sets
from hahn_pool import fit
KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'Darkwave Disco'
SHOTS = 1200; SEED = 20261023; SETS_USED = (0,)
TAUS = {'T1': [0, 20, 50, 100, 170, 260, 400, 600], 'RAMSEY': [0, 5, 10, 20, 35, 60, 100, 160],
        'HAHN': [0, 20, 40, 90, 160, 220, 300, 400], 'CPMG4': [0, 40, 90, 160, 260, 400, 600, 800], 'CPMG16': [0, 40, 90, 160, 260, 400, 600, 800]}
NPULSE = {'T1': 1, 'RAMSEY': 0, 'HAHN': 1, 'CPMG4': 4, 'CPMG16': 16}
ECHO = ('HAHN', 'CPMG4', 'CPMG16')

def circuit(n, S, rung, t, quad=0):
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(n, len(S))
    if rung == 'T1':
        qc.x(S); qc.barrier()
        if t: qc.delay(t, S, unit='us')
    elif rung == 'RAMSEY':
        qc.sx(S); qc.barrier()
        if t: qc.delay(t, S, unit='us')
        qc.barrier()
        if quad: qc.rz(np.pi / 2, S)
        qc.sxdg(S)
    else:
        N = {'HAHN': 1, 'CPMG4': 4, 'CPMG16': 16}[rung]; g = qc.x if rung == 'HAHN' else qc.y
        qc.sx(S); qc.barrier()
        for k in range(N):
            w = t / (2 * N) if k == 0 else t / N
            if w: qc.delay(w, S, unit='us')
            qc.barrier(); g(S); qc.barrier()
        if t: qc.delay(t / (2 * N), S, unit='us')
        qc.barrier(); qc.sxdg(S)
    qc.barrier(); qc.measure(S, list(range(len(S))))
    return qc

def items(n, sets):
    out = []
    for si in SETS_USED:
        S = sets[si]
        for rung, ts in TAUS.items():
            for t in ts:
                for quad in ((0, 1) if rung == 'RAMSEY' else (0,)):
                    for cp in (0, 1): out.append((circuit(n, S, rung, t, quad), dict(set=si, rung=rung, tau_us=t, quad=quad, copy=cp)))
    random.Random(SEED).shuffle(out); return out

def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)

# ---------- analysis ----------
def curves(meta, bits_of, sets):
    """per qubit, per rung: copy -> {tau: signal}. Signal: T1 P1; HAHN P1; CPMG 1-P1; RAMSEY 0.5 + c/2 with the shot-noise floor removed."""
    raw = {}
    for i, m in enumerate(meta):
        S = sets[m['set']]; b = bits_of(i)
        for j, q in enumerate(S): raw.setdefault(q, {}).setdefault(m['rung'], {}).setdefault(m['copy'], {})[(m['tau_us'], m['quad'])] = (float(b[:, j].mean()), b.shape[0])
    out = {}
    for q, R in raw.items():
        for rung, cps in R.items():
            for cp, d in list(cps.items()) + [('pool', None)]:
                if cp == 'pool':
                    keys = d0 = cps[0].keys()
                    d = {k: ((cps[0][k][0] + cps[1][k][0]) / 2, cps[0][k][1] + cps[1][k][1]) for k in keys}
                sig = {}
                for t in TAUS[rung]:
                    if rung == 'RAMSEY':
                        (pa, na), (pb, nb) = d[(t, 0)], d[(t, 1)]; xa, xb = 1 - 2 * pa, 1 - 2 * pb
                        c2 = xa ** 2 + xb ** 2 - (1 - xa ** 2) / na - (1 - xb ** 2) / nb; sig[t] = 0.5 + np.sqrt(max(c2, 0)) / 2
                    else:
                        p = d[(t, 0)][0]; sig[t] = 1 - p if rung in ('CPMG4', 'CPMG16') else p
                out.setdefault(q, {}).setdefault(rung, {})[cp] = sig
    return out

def times(cv, free):
    """per qubit, per rung: (T, A) from the pooled copies, and per copy. T1 always free-base (decays to the readout floor, not 1/2)."""
    T = {}
    for q, R in cv.items():
        for rung, cps in R.items():
            for cp, sig in cps.items():
                ts = TAUS[rung]; y = [sig[t] for t in ts]
                if rung == 'T1': tt, A, _ = fit(ts, y, True); ok = A > 0.4
                else: tt, A, _ = fit(ts, y, free if rung != 'RAMSEY' else False); ok = A > 0.2
                T.setdefault(q, {}).setdefault(rung, {})[cp] = tt if ok and np.isfinite(tt) and tt < 1990 else np.nan
    return T

def med_ci(x, rng, B=4000):
    x = np.asarray([v for v in x if np.isfinite(v)])
    if len(x) < 10: return np.nan, np.nan, np.nan, len(x)
    m = float(np.median(x)); lo, hi = np.percentile([np.median(rng.choice(x, len(x))) for _ in range(B)], [2.5, 97.5]); return m, lo, hi, len(x)
def four(lo, hi, ref=1.0, gap_shift=False):
    if gap_shift or not np.isfinite(lo): return 'CANNOT_TELL'
    return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'
def copygap(T, rung, rng):
    g = [T[q][rung][0] / T[q][rung][1] for q in T if rung in T[q] and np.isfinite(T[q][rung].get(0, np.nan)) and np.isfinite(T[q][rung].get(1, np.nan))]
    m, lo, hi, n = med_ci(g, rng); return m, lo, hi, (lo > 1.05 or hi < 0.95)   # shifted = clears 1 by more than 5 % (synthetic floor 3.2 %)

def readings(cv, stated, label):
    rng = np.random.default_rng(1); res = {}
    for free in (False, True):
        T = times(cv, free); P = {q: {r: T[q][r]['pool'] for r in T[q]} for q in T}; fl = 'FREE ' if free else 'FIXED'
        print(f'-- {label} | {fl} fit (T1 always free-base; Ramsey on the two-quadrature magnitude)')
        shifted = {r: copygap(T, r, rng) for r in TAUS}
        for r in TAUS:
            m, lo, hi, sh = shifted[r]; print(f'   copy gap {r:7} T(copy0)/T(copy1) median {m:.3f} [{lo:.3f}, {hi:.3f}]' + ('  SHIFTED -> rung CANNOT_TELL' if sh else ''))
        med = {r: np.nanmedian([P[q][r] for q in P if r in P[q]]) for r in TAUS}
        print('   chip median times (us): ' + ', '.join(f'{r} {med[r]:.0f}' for r in TAUS))
        # R1, R2: delivered / IBM stated, same record
        for r, key in (('T1', 'T1'), ('HAHN', 'T2')):
            m, lo, hi, n = med_ci([P[q][r] / stated[key][q] for q in P if q in stated[key] and r in P[q]], rng)
            out = four(lo, hi, 1, shifted[r][3]); res[(fl, 'R1' if r == 'T1' else 'R2')] = out
            print(f'   {"R1" if r == "T1" else "R2"} delivered/stated {key:3} median {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {out}')
        # R3: the ladder, per qubit ratio to the Hahn rung
        for r in ('RAMSEY', 'CPMG4', 'CPMG16'):
            m, lo, hi, n = med_ci([P[q][r] / P[q]['HAHN'] for q in P if r in P[q] and 'HAHN' in P[q]], rng)
            out = four(lo, hi, 1, shifted[r][3] or shifted['HAHN'][3]); res[(fl, 'R3', r)] = out
            print(f'   R3 T2({r})/T2(HAHN) per qubit median {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {out}')
        # R4: distance from the T1 limit, T2/(2 T1) per qubit, every echo rung; KNOWN-GOOD bound on HAHN
        for r in ECHO:
            m, lo, hi, n = med_ci([P[q][r] / (2 * P[q]['T1']) for q in P if r in P[q] and 'T1' in P[q]], rng)
            out = four(lo, hi, 1); res[(fl, 'R4', r)] = (m, lo, hi, out)
            print(f'   R4 T2({r})/(2 T1) per qubit median {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {out}')
        # NOISE: the class, from R3 CPMG16 and R4 (rule in the plan)
        c16 = res[(fl, 'R3', 'CPMG16')]; h4 = res[(fl, 'R4', 'HAHN')]; c4 = res[(fl, 'R4', 'CPMG16')]
        if c16 == 'CANNOT_TELL' or h4[3] == 'CANNOT_TELL': cls = 'CANNOT_TELL'
        elif h4[2] >= 0.9: cls = 'T1-LIMITED (Hahn already within 10 % of 2 T1)'
        elif c16 == 'ABOVE' and c4[3] == 'BELOW': cls = 'SLOW NOISE (decoupling lengthens T2, still below 2 T1)'
        elif c16 == 'ABOVE': cls = 'SLOW NOISE, REMOVED (decoupling reaches 2 T1)'
        elif c16 == 'CONSISTENT': cls = 'WHITE (decoupling does not lengthen T2; Hahn below 2 T1)'
        else: cls = 'BELOW: decoupling SHORTENS T2 (pulse error or fast noise)'
        res[(fl, 'NOISE')] = cls; print(f'   NOISE class: {cls}')
    return res

def synth(world, seed=5, nq=90):
    """Planted worlds. T1 per qubit 120-320 us; stated = true except 'short'. Echo T2 per rung from the world:
       slow    : Tphi(N) grows as N^0.4 (1/f-like): Hahn below 2T1, CPMG longer
       white   : Tphi the same for every N
       t1limit : no dephasing, T2 = 2 T1 on every echo rung
       mixed   : half the qubits slow, half white
       short   : as slow, but every true T1 and T2 = 0.85 x stated
       pulse   : as white, but each pi pulse loses 0.4 % contrast (an amplitude loss: must NOT read as a T2 change)
       drift   : as slow, but copy 1 of every echo rung has T2 x 0.80 (the copies disagree: echo rungs must read CANNOT_TELL)"""
    rng = np.random.default_rng(seed); qs = list(range(nq)); T1 = {q: float(rng.uniform(120, 320)) for q in qs}
    stated = {'T1': dict(T1), 'T2': {}}; tp0 = {q: float(rng.uniform(150, 600)) for q in qs}
    meta, bits = [], []
    sets = [qs]; f = 0.85 if world == 'short' else 1.0
    for q in qs: stated['T2'][q] = 1 / (1 / (2 * T1[q]) + 1 / tp0[q])
    def tphi(q, N):
        if world == 't1limit': return np.inf
        kind = world if world != 'mixed' else ('slow' if q % 2 else 'white')
        if kind == 'drift': kind = 'slow'
        if kind in ('slow', 'short'): return tp0[q] * N ** 0.4
        return tp0[q]
    for rung, ts in TAUS.items():
        for t in ts:
            for quad in ((0, 1) if rung == 'RAMSEY' else (0,)):
                for cp in (0, 1):
                    p = []
                    for q in qs:
                        t1 = T1[q] * f
                        if rung == 'T1': pr = 0.02 + 0.95 * np.exp(-t / t1)
                        elif rung == 'RAMSEY':
                            tst = 1 / (1 / (2 * t1) + 1 / (0.4 * tphi(q, 1) * f)) if np.isfinite(tphi(q, 1)) else 2 * t1
                            phi = 2 * np.pi * 0.015 * t * (1 + q % 3)
                            c = 0.95 * np.exp(-t / tst); pr = 0.5 - 0.5 * c * (np.cos(phi) if quad == 0 else np.sin(phi))
                        else:
                            N = NPULSE[rung]; tp = tphi(q, N) * f; T2 = 1 / (1 / (2 * t1) + (1 / tp if np.isfinite(tp) else 0))
                            if world == 'drift' and cp == 1: T2 *= 0.80
                            loss = (1 - 0.004) ** N if world == 'pulse' else 1.0; A = 0.47 * loss * np.exp(-t / T2)
                            pr = 0.5 + A if rung == 'HAHN' else 0.5 - A
                        p.append(np.clip(pr, 0, 1))
                    meta.append(dict(set=0, rung=rung, tau_us=t, quad=quad, copy=cp)); bits.append(rng.random((SHOTS, nq)) < np.array(p))
    return meta, bits, sets, stated

def stated_of(props):
    st = {'T1': {}, 'T2': {}}
    for q, ps in enumerate(props['qubits']):
        for x in ps:
            if x['name'] in st and x['value'] > 0: st[x['name']][q] = x['value']
    return st

if __name__ == '__main__':
    if '--synth' in sys.argv:
        w = sys.argv[sys.argv.index('--synth') + 1]; meta, bits, sets, stated = synth(w)
        print(f'=== synthetic world {w}: {synth.__doc__.split(w + " ")[1].split(chr(10))[0].strip(": ") if (w + " ") in synth.__doc__ else ""}')
        readings(curves(meta, lambda i: bits[i], sets), stated, f'synth {w}'); sys.exit(0)
    if '--fetch' in sys.argv:
        chip = sys.argv[sys.argv.index('--fetch') + 1]; rec = json.load(open(f'{HERE}/job_ladder_{chip}.json'))
        res = service().job(rec['job_id']).result(); out = {f'c_{i}': pr.data.c.to_bool_array(order='little') for i, pr in enumerate(res)}
        np.savez_compressed(f'{HERE}/raw_ladder_{chip}.npz', **out); print('saved', chip, len(res), 'pubs'); sys.exit(0)
    if '--analyze' in sys.argv:
        chip = sys.argv[sys.argv.index('--analyze') + 1]
        rec = json.load(open(f'{HERE}/job_ladder_{chip}.json')); z = np.load(f'{HERE}/raw_ladder_{chip}.npz'); props = json.load(open(f'{HERE}/props_ladder_{chip}.json'))
        stated = stated_of(props); cv = curves(rec['meta'], lambda i: z[f'c_{i}'][:, :len(rec['sets'][rec['meta'][i]['set']])], rec['sets'])
        h0 = np.median([cv[q]['HAHN']['pool'][0] for q in cv]); r0 = np.median([2 * cv[q]['RAMSEY']['pool'][0] - 1 for q in cv])
        print(f'=== {chip} job {rec["job_id"]}; IBM record {props["last_update_date"]}; submitted {rec["submitted_utc"]}')
        print(f'KNOWN-GOOD (1) median Hahn P(1) at t = 0: {h0:.3f} (>= 0.90): {"MET" if h0 >= 0.90 else "NOT MET"} | (1b) median Ramsey coherence at t = 0: {r0:.3f} (>= 0.85): {"MET" if r0 >= 0.85 else "NOT MET"}')
        res = readings(cv, stated, chip)
        for fl in ('FIXED', 'FREE '):
            m, lo, hi, out = res[(fl, 'R4', 'HAHN')]
            print(f'KNOWN-GOOD (2) {fl}: T2(HAHN)/(2 T1) upper bound {hi:.3f} must be <= 1.05 (no qubit population may beat the T1 limit): {"MET" if hi <= 1.05 else "NOT MET"}')
        sys.exit(0)
    CHIP = sys.argv[1]; FIRE = '--fire' in sys.argv
    from qiskit import transpile
    s = service(); b = s.backend(CHIP); n = b.num_qubits
    edges = sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}); sets = sparse_sets(n, edges); its = items(n, sets)
    tq = transpile([c for c, _ in its], b, optimization_level=1, initial_layout=list(range(n)), scheduling_method='alap', seed_transpiler=7)
    bad = 0
    for c, (_, m) in zip(tq, its):     # every pi pulse survived transpilation (a cancelled pair would fake a longer T2)
        ops = c.count_ops(); S = sets[m['set']]; k = NPULSE[m['rung']]
        npi = ops.get('x', 0) + ops.get('y', 0)
        if m['rung'] in ('T1', 'HAHN', 'CPMG4', 'CPMG16') and npi < k * len(S): bad += 1
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: set sizes {[len(sets[i]) for i in SETS_USED]}, {len(tq)} circuits, shots {SHOTS}, pulse-count failures {bad}, ESTIMATED {est:.0f} s, remaining {rem} ({INSTANCE})')
    print('   ops in one CPMG16 circuit:', dict(tq[[i for i, (_, m) in enumerate(its) if m['rung'] == 'CPMG16' and m['tau_us'] == 0][0]].count_ops()))
    if not FIRE: sys.exit(0)
    if bad: print('STOP: pulses lost in transpilation'); sys.exit(1)
    if rem is not None and rem < est + 10: print('STOP: not enough seconds'); sys.exit(1)
    from qiskit_ibm_runtime import SamplerV2
    json.dump(b.properties().to_dict(), open(f'{HERE}/props_ladder_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), meta=[m for _, m in its], sets=sets, shots=SHOTS, keyfile=KEYFILE, instance=INSTANCE)
    json.dump(rec, open(f'{HERE}/job_ladder_{CHIP}.json', 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
