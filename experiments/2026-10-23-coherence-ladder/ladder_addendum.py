"""ladder_addendum.py -- ADDENDUM_LADDER.md, written after the jobs were fetched and BEFORE any result was read (/q pass, 23 Oct 2026).
ladder.py is unchanged (its hash holds); every registered reading still prints from ladder.py --analyze. This file adds:
  --post CHIP     save IBM's record as it stands now (props_ladder_CHIP_post.json) and the job's run times (zero QPU seconds)
  --synth WORLD [SEED]   extra planted worlds: gauss (Gaussian echo decay whose T grows with N), white/slow seeds for the floor
  --floor         false-resolution floor: 10 seeds each of white and slow
  --analyze CHIP  the addendum readings (below)
Addendum readings (both fits):
  X1  R3b = T2(CPMG16)/T2(CPMG4): both rungs use Y pulses, so this ratio changes ONLY the pulse count (R3 against HAHN also changes the
      pulse axis X -> Y). Rule: the NOISE class says SLOW only if R3b is not BELOW; otherwise the class is CANNOT_TELL (axis-confounded).
  X2  Ramsey refit on waits where the coherence is at least 5 x its shot-noise floor; if its R3 outcome differs from the registered one,
      RAMSEY is CANNOT_TELL.
  X3  per rung: qubits dropped (amplitude cut) and capped (T >= 1990 us), and the CPMG known-good: median 1-P(1) at wait 0 >= 0.90.
  X4  T1 shape: per qubit, the 1/e crossing of the pooled T1 curve (linear interpolation, floor removed) over the fitted T1; IBM's docs
      call T1 one time constant (claims row 24). Descriptive: median and interval, no bill reading.
  X5  record drift: R1 and R2 recomputed against the record as it stands after the job; count of qubits whose T1/T2 value changed.
  X6  schedule: per circuit, transpiled duration minus nominal wait, by rung; max spread inside a rung (pulse time is not idle time).
  TWO-CHIP RULE: an R3, R3b or NOISE outcome is a finding only if both chips read it; one chip is printed as one chip."""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ladder as L

def synth2(world, seed):
    rng = np.random.default_rng(seed); nq = 90; qs = list(range(nq)); T1 = {q: float(rng.uniform(120, 320)) for q in qs}
    tp0 = {q: float(rng.uniform(150, 600)) for q in qs}; stated = {'T1': dict(T1), 'T2': {q: 1 / (1 / (2 * T1[q]) + 1 / tp0[q]) for q in qs}}
    meta, bits = [], []
    for rung, ts in L.TAUS.items():
        for t in ts:
            for quad in ((0, 1) if rung == 'RAMSEY' else (0,)):
                for cp in (0, 1):
                    p = []
                    for q in qs:
                        if rung == 'T1': pr = 0.02 + 0.95 * np.exp(-t / T1[q])
                        elif rung == 'RAMSEY':
                            tst = 1 / (1 / (2 * T1[q]) + 1 / (0.4 * tp0[q])); phi = 2 * np.pi * 0.015 * t * (1 + q % 3); c = 0.95 * np.exp(-t / tst)
                            pr = 0.5 - 0.5 * c * (np.cos(phi) if quad == 0 else np.sin(phi))
                        else:
                            N = L.NPULSE[rung]; tp = tp0[q] * (N ** 0.4 if world in ('slow', 'gauss') else 1)
                            dec = np.exp(-t / (2 * T1[q])) * (np.exp(-(t / tp) ** 2) if world == 'gauss' else np.exp(-t / tp))
                            A = 0.47 * dec; pr = 0.5 + A if rung == 'HAHN' else 0.5 - A
                        p.append(np.clip(pr, 0, 1))
                    meta.append(dict(set=0, rung=rung, tau_us=t, quad=quad, copy=cp)); bits.append(rng.random((L.SHOTS, nq)) < np.array(p))
    return meta, bits, [qs], stated

def extra(cv, stated, post=None, label=''):
    rng = np.random.default_rng(2); out = {}
    for free in (False, True):
        fl = 'FREE ' if free else 'FIXED'; T = L.times(cv, free); P = {q: {r: T[q][r]['pool'] for r in T[q]} for q in T}
        m, lo, hi, n = L.med_ci([P[q]['CPMG16'] / P[q]['CPMG4'] for q in P], rng); x1 = L.four(lo, hi); out[(fl, 'X1')] = x1
        print(f'   X1 {fl} T2(CPMG16)/T2(CPMG4) median {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {x1}')
        # X2 Ramsey censored refit
        r2 = []
        for q in cv:
            sig = cv[q]['RAMSEY']['pool']; ts = [t for t in L.TAUS['RAMSEY'] if (2 * sig[t] - 1) >= 5 * np.sqrt(2 / (2 * L.SHOTS))]
            if len(ts) >= 4:
                tt, A, _ = L.fit(ts, [sig[t] for t in ts], False)
                if A > 0.2 and tt < 1990 and np.isfinite(P[q].get('HAHN', np.nan)): r2.append(tt / P[q]['HAHN'])
        m, lo, hi, n = L.med_ci(r2, rng); print(f'   X2 {fl} Ramsey (waits >= 5x floor) / HAHN median {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {L.four(lo, hi)}')
        out[(fl, 'X2')] = L.four(lo, hi)
        # X3 dropped / capped
        print('   X3 ' + fl + ' dropped/capped per rung: ' + ', '.join(f'{r} {sum(1 for q in T if not np.isfinite(T[q][r]["pool"]))}' for r in L.TAUS))
        if post is not None:
            for r, key in (('T1', 'T1'), ('HAHN', 'T2')):
                m, lo, hi, n = L.med_ci([P[q][r] / post[key][q] for q in P if q in post[key]], rng)
                print(f'   X5 {fl} delivered/stated {key} against the POST-job record: {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {L.four(lo, hi)}')
    for r in ('CPMG4', 'CPMG16'):
        v = np.median([cv[q][r]['pool'][0] for q in cv]); print(f'   X3 known-good {r} median 1-P(1) at wait 0: {v:.3f} (>= 0.90): {"MET" if v >= 0.90 else "NOT MET"}')
    # X4 T1 shape
    rat = []
    for q in cv:
        sig = cv[q]['T1']['pool']; ts = np.array(L.TAUS['T1'], float); y = np.array([sig[t] for t in L.TAUS['T1']])
        tt, A, _ = L.fit(ts, y, True); B = y[0] - A
        if A > 0.4 and np.isfinite(tt):
            target = B + A / np.e; k = np.where(y <= target)[0]
            if len(k) and k[0] > 0:
                i = k[0]; te = ts[i - 1] + (y[i - 1] - target) * (ts[i] - ts[i - 1]) / (y[i - 1] - y[i]); rat.append(te / tt)
    m, lo, hi, n = L.med_ci(rat, rng); print(f'   X4 T1 shape: 1/e crossing / fitted T1 median {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {L.four(lo, hi)} (1 = one time constant)')
    return out

if __name__ == '__main__':
    a = sys.argv
    if '--synth' in a:
        w = a[a.index('--synth') + 1]; seed = int(a[a.index('--synth') + 2]) if len(a) > a.index('--synth') + 2 else 5
        meta, bits, sets, st = synth2(w, seed); cv = L.curves(meta, lambda i: bits[i], sets)
        print(f'=== addendum synthetic {w} seed {seed}'); res = L.readings(cv, st, f'{w} {seed}'); extra(cv, st, None, w); sys.exit(0)
    if '--post' in a:
        chip = a[a.index('--post') + 1]; s = L.service(); rec = json.load(open(f'{HERE}/job_ladder_{chip}.json'))
        b = s.backend(chip); json.dump(b.properties().to_dict(), open(f'{HERE}/props_ladder_{chip}_post.json', 'w'), default=str)
        j = s.job(rec['job_id']); mt = j.metrics(); ts = mt.get('timestamps', {}); u = mt.get('usage', {})
        json.dump({'timestamps': ts, 'usage': u}, open(f'{HERE}/run_ladder_{chip}.json', 'w'), default=str)
        print(chip, 'submitted', rec['submitted_utc'], '| timestamps', ts, '| usage', u); sys.exit(0)
    if '--analyze' in a:
        chip = a[a.index('--analyze') + 1]; rec = json.load(open(f'{HERE}/job_ladder_{chip}.json')); z = np.load(f'{HERE}/raw_ladder_{chip}.npz')
        st = L.stated_of(json.load(open(f'{HERE}/props_ladder_{chip}.json'))); post = L.stated_of(json.load(open(f'{HERE}/props_ladder_{chip}_post.json')))
        ch = {k: sum(1 for q in st[k] if q in post[k] and post[k][q] != st[k][q]) for k in st}
        print(f'=== ADDENDUM {chip} | X5 qubits whose IBM value changed between submission and the post-job record: T1 {ch["T1"]}, T2 {ch["T2"]}')
        cv = L.curves(rec['meta'], lambda i: z[f'c_{i}'][:, :len(rec['sets'][rec['meta'][i]['set']])], rec['sets']); extra(cv, st, post, chip); sys.exit(0)

def sched(chip):
    from qiskit import transpile
    s = L.service(); b = s.backend(chip); n = b.num_qubits
    edges = sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}); sets = L.sparse_sets(n, edges); its = L.items(n, sets)
    tq = transpile([c for c, _ in its], b, optimization_level=1, initial_layout=list(range(n)), scheduling_method='alap', seed_transpiler=7)
    by = {}
    for c, (_, m) in zip(tq, its): by.setdefault(m['rung'], []).append((m['tau_us'], c.estimate_duration(b.target, 's') * 1e6 - m['tau_us']))
    for r, v in by.items():
        ex = np.array([x for _, x in v]); print(f'   X6 {chip} {r:7} duration minus nominal wait: {ex.min():.2f} to {ex.max():.2f} us (spread {ex.max() - ex.min():.2f})')
if __name__ == '__main__' and '--sched' in sys.argv: sched(sys.argv[sys.argv.index('--sched') + 1])
