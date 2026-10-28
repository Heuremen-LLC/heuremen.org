#!/usr/bin/env python3
"""verify_cost.py -- recompute the cost of a mid-circuit measurement from the reduced data.

    python3 verify_cost.py

Standard library only. No account, no network. Reads two files in data/:
  fez_measurement_series_20261011.json   three IBM jobs on ibm_fez, 11 Oct 2026 (1 ... 128 measurements in a row)
  kingston_spectators_20261004_06.json   six IBM jobs on ibm_kingston, 4 and 6 Oct 2026 (1 ... 32 measurements)
and prints every recomputed figure with a plain label.
"""
import json, math, os, random, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
FEZ = json.load(open(os.path.join(HERE, 'data', 'fez_measurement_series_20261011.json')))
KING = json.load(open(os.path.join(HERE, 'data', 'kingston_spectators_20261004_06.json')))

# field names exactly as stored in the data files
F_LADDER = 'readout_counts'             # numbers of consecutive measurements in the ladder
F_BITS = 'measured_qubit_bits'            # per-qubit records for the measured qubits (ibm_fez jobs)
F_ON_K = 'level_after_k_readouts'          # measured-qubit level after k measurements
F_SPEC_K = 'qubit_level_after_k_readouts'      # spectator level after k measurements
F_STORED = 'stored_fit' # per-fit bootstrap summary stored with each ibm_fez job
F_STATED_US = 'stated_readout_us'       # stated length of one measurement
F_MEAS_US = 'measured_readout_us'  # measured length of one measurement, per day


def report(label, got, note=''):
    go = f'{got:.4g}' if isinstance(got, float) else str(got)
    print(f'  {label:<72} {go}' + (f'   [{note}]' if note else ''))


# ---------------------------------------------------------------- small numerical tools

def median(v): return st.median(v)
def mean(v): return sum(v) / len(v)
def sd(v): return st.stdev(v)


def interp(x, xs, ys):
    """piecewise-linear, clamped"""
    if x <= xs[0]: return ys[0]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            f = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + f * (ys[i] - ys[i - 1])
    return ys[-1]


def loginterp(x, xs, ys): return math.exp(interp(x, xs, [math.log(y) for y in ys]))


def stretched_fit(x, y, s):
    """best y = A exp(-(x/k0)^b); A solved in closed form at each (k0, b). Grid search:
    coarse grid b in 0.2 .. 3.0 (141 steps) x k0 log-spaced 0.5 .. 4000 (241 steps), then five local refinements."""
    w = [1.0 / v ** 2 for v in s]

    def chi(k0, b):
        p = [math.exp(-(xi / k0) ** b) for xi in x]
        den = sum(wi * pi * pi for wi, pi in zip(w, p))
        if den <= 0: return float('inf'), 0.0
        A = sum(wi * pi * yi for wi, pi, yi in zip(w, p, y)) / den
        return sum(wi * (yi - A * pi) ** 2 for wi, pi, yi in zip(w, p, y)), A

    best = None
    lo, hi = math.log(0.5), math.log(4000.0)
    for ib in range(141):
        b = 0.2 + 2.8 * ib / 140
        for ik in range(241):
            k0 = math.exp(lo + (hi - lo) * ik / 240)
            c, A = chi(k0, b)
            if best is None or c < best[0]: best = (c, k0, b, A)
    c, k0, b, A = best
    for _ in range(5):
        k0_round = k0                                    # the k0 grid is fixed at the start of each round
        for i in range(41):                              # the b grid is re-centred on the current b at each step of the outer loop
            k0t = k0_round * (0.85 + 0.33 * i / 40)
            b_now = b
            for j in range(41):
                bt = b_now * (0.85 + 0.33 * j / 40)
                ct, At = chi(k0t, bt)
                if ct < c: c, k0, b, A = ct, k0t, bt, At
    return A, k0, b, c


def local_fit(x, y, s, start):
    """cheap re-fit for the resampling check: pattern search around a starting (k0, b)."""
    w = [1.0 / v ** 2 for v in s]

    def chi(k0, b):
        if k0 <= 0.1 or b <= 0.1 or b > 4: return float('inf'), 0.0
        p = [math.exp(-(xi / k0) ** b) for xi in x]
        den = sum(wi * pi * pi for wi, pi in zip(w, p))
        A = sum(wi * pi * yi for wi, pi, yi in zip(w, p, y)) / den
        return sum(wi * (yi - A * pi) ** 2 for wi, pi, yi in zip(w, p, y)), A

    k0, b = start; c, _ = chi(k0, b); step = 0.2
    while step > 1e-5:
        moved = False
        for dk, db in ((1 + step, 1), (1 / (1 + step), 1), (1, 1 + step), (1, 1 / (1 + step))):
            ct, _ = chi(k0 * dk, b * db)
            if ct < c: c, k0, b, moved = ct, k0 * dk, b * db, True
        if not moved: step /= 2
    return k0, b


def runs_p(signs):
    """exact Wald-Wolfowitz P(runs <= observed)"""
    s = [v for v in signs if v != 0]; n = len(s); n1 = sum(1 for v in s if v > 0); n2 = n - n1
    if n1 == 0 or n2 == 0: return 1.0
    obs = 1 + sum(1 for i in range(1, n) if s[i] != s[i - 1])

    def arr(r):
        if r % 2 == 0:
            k = r // 2; return 2 * math.comb(n1 - 1, k - 1) * math.comb(n2 - 1, k - 1)
        k = (r - 1) // 2
        return math.comb(n1 - 1, k) * math.comb(n2 - 1, k - 1) + math.comb(n1 - 1, k - 1) * math.comb(n2 - 1, k)
    return min(1.0, sum(arr(r) for r in range(2, obs + 1)) / math.comb(n, n1))


# ================================================================ PART B: the 11 October series on ibm_fez

def fez_job(j):
    """the two curves of one job: the measured qubit's, and the spectators'."""
    ks = FEZ[F_LADDER]
    bits = j[F_BITS]; specs = j['spectators']
    on = [mean([b[F_ON_K][i] for b in bits]) for i in range(len(ks))]
    plain = [mean([b['on_after_plain_wait'][i] for b in bits]) for i in range(len(FEZ['own_waits_us']))]
    hand = [median([c[F_SPEC_K][i] for c in specs]) for i in range(len(ks))]
    laid = [median([c['hand_after_laid_wait'][i] for c in specs]) for i in range(len(FEZ['laid_waits_us']))]
    return on, plain, hand, laid, len(bits), len(specs)


def seat_fit(raw, floor, widths, ks):
    span = raw[0] - floor
    r = [(v - floor) / span for v in raw]
    # first point is the normaliser (width 0) and is excluded from the fit
    x = ks[1:]; y = r[1:]; s = widths[1:]
    A, k0, b, chi2 = stretched_fit(x, y, s)
    fit = [A * math.exp(-(xi / k0) ** b) for xi in x]
    dof = max(1, len(x) - 3)
    runs = runs_p([1 if yi - fi > 0 else -1 if yi - fi < 0 else 0 for yi, fi in zip(y, fit)])
    return dict(A=A, k0=k0, b=b, chi2dof=chi2 / dof, runs_p=runs, span=span, r=r, x=x, y=y, s=s)


def part_b():
    print('\nB.  ibm_fez, 11 October, three jobs: ' + ', '.join(f"run {j['run']} = IBM job {j['job_id']}" for j in FEZ['jobs']))
    ks = FEZ[F_LADDER]; res = []
    for j in FEZ['jobs']:
        on, plain, hand, laid, nb, nc = fez_job(j)
        pr = j[F_STORED]
        s1 = seat_fit(on, 0.0, [0.0] + pr['seat1_measured_qubit']['bootstrap_width_of_normalised_level'][1:], ks)
        s2 = seat_fit(hand, laid[-1], [0.0] + pr['seat2_spectators']['bootstrap_width_of_normalised_level'][1:], ks)
        res.append(dict(job=j, on=on, plain=plain, hand=hand, laid=laid, nb=nb, nc=nc, s1=s1, s2=s2))
    print('  (curve fits are re-run here from the stored per-qubit levels; the shot-level widths are read from the file)')

    # --- shape exponents
    b1 = [r['s1']['b'] for r in res]; b2 = [r['s2']['b'] for r in res]
    print('\n  shape exponent b per job  (measured qubit | spectators):')
    for r in res:
        pr = r['job'][F_STORED]
        print(f"      run {r['job']['run']}:  {r['s1']['b']:.3f} | {r['s2']['b']:.3f}     "
              f"(stored fit values {pr['seat1_measured_qubit']['stored_b']:.3f} | {pr['seat2_spectators']['stored_b']:.3f};"
              f" per-fit bootstrap widths +-{pr['seat1_measured_qubit']['b_width']:.3f} | +-{pr['seat2_spectators']['b_width']:.3f})")
    m1, w1, m2, w2 = mean(b1), sd(b1), mean(b2), sd(b2)
    report('measured-qubit exponent, mean of 3 jobs', m1)
    report('  its width (sample SD of 3 jobs)', w1, 'n = 3, a spread of three numbers')
    report('spectator exponent, mean of 3 jobs', m2)
    report('  its width (sample SD of 3 jobs)', w2, 'n = 3')
    report('separation of the two, in combined widths', (m2 - m1) / math.hypot(w1, w2), 'divides by the spread of single jobs')
    report('measured qubit below a plain exponential, widths', (1 - m1) / w1)
    report('spectators above a plain exponential, widths', (m2 - 1) / w2)
    sem1, sem2 = w1 / math.sqrt(3), w2 / math.sqrt(3)
    report('same separation using the uncertainty OF THE MEAN (SD/sqrt3)', (m2 - m1) / math.hypot(sem1, sem2))

    # --- chi2 / dof
    c1 = [r['s1']['chi2dof'] for r in res]; c2 = [r['s2']['chi2dof'] for r in res]
    print('\n  chi2/dof per job: measured qubit ' + ' / '.join(f'{v:.2f}' for v in c1) + '   spectators ' + ' / '.join(f'{v:.2f}' for v in c2))
    report('measured qubit chi2/dof (three jobs)', ' / '.join(f'{v:.2f}' for v in c1))
    report('spectators chi2/dof (three jobs)', ' / '.join(f'{v:.2f}' for v in c2))
    report('worst of the six chi2/dof', max(c1 + c2))
    report('run-structure p, measured qubit (three jobs)', '/'.join(f"{r['s1']['runs_p']:.2f}" for r in res))
    report('run-structure p, spectators (three jobs)', '/'.join(f"{r['s2']['runs_p']:.2f}" for r in res))
    c3 = [r['job'][F_STORED]['seat3_no_measurement_control']['stored_chi2_per_dof'] for r in res]
    print('  (the third arm of each job, a no-measurement control with 5 points, is not one of the six fits and does not fit:'
          ' stored chi2/dof ' + ' / '.join(f'{v:.1f}' for v in c3) + ')')
    # --- resampling check of the width of b from the stored point widths
    random.seed(20261011)
    print('\n  cross-check of the single-fit widths of b by redrawing each point from its stored width (60 draws):')
    for r in res:
        for nm in ('s1', 's2'):
            f = r[nm]; bs = []
            for _ in range(60):
                yy = [yi + random.gauss(0, si) for yi, si in zip(f['y'], f['s'])]
                bs.append(local_fit(f['x'], yy, f['s'], (f['k0'], f['b']))[1])
            pw = r['job'][F_STORED]['seat1_measured_qubit' if nm == 's1' else 'seat2_spectators']['b_width']
            print(f"      run {r['job']['run']} {'measured qubit' if nm == 's1' else 'spectators    '}: spread {sd(bs):.3f}   (shot-level bootstrap {pw:.3f})")

    # --- headline curve numbers (job 1)
    r1 = res[0]
    print('\n  measured-qubit levels (job 1, run %d):' % r1['job']['run'])
    report('level on the F0-F1 axis after 1 measurement', r1['on'][0])
    report('level after 128 measurements', r1['on'][-1])
    report('level after 40 us of plain idle (same job, same bits)', r1['plain'][-1])
    cross = []
    for r in res:
        o = r['on']
        for i in range(1, len(ks)):
            if o[i] <= 0.5 < o[i - 1]:
                cross.append(ks[i - 1] + (ks[i] - ks[i - 1]) * (o[i - 1] - 0.5) / (o[i - 1] - o[i])); break
    report('number of measurements at which level = 1/2 (job 1)', cross[0], 'straight line between the 24 and 32 points')
    print(f"      the same crossing in jobs 2 and 3: {cross[1]:.1f} and {cross[2]:.1f}")
    print('      level after 1 / 128 measurements / 40 us idle, all three jobs:')
    for r in res:
        print(f"        run {r['job']['run']}: {r['on'][0]:.4f} / {r['on'][-1]:.4f} / {r['plain'][-1]:.4f}")

    # --- is the 40 us idle an equal-time idle?
    L = FEZ[F_STATED_US]
    print(f'\n  EQUAL-TIME COMPARISON?  128 measurements x {L} us stated length = {128 * L:.0f} us, not 40 us.')
    report('time taken by 128 measurements (stated length), us', 128 * L, 'the 40 us idle arm is 5x shorter')
    idle_t = FEZ['own_waits_us']; pairs = []
    print('  equal-time pairs in the same job (measured qubit; k measurements vs an idle of k x stated length):')
    for k, t in ((6, 10), (12, 20), (24, 40)):
        i = ks.index(k); w = idle_t.index(t)
        row = [(r['on'][i], r['plain'][w]) for r in res]
        pairs.append((k, t, row))
        print(f"      k={k:>2} ({k * L:5.1f} us) vs idle {t} us:  " + '   '.join(f'{a:.3f} vs {b:.3f}' for a, b in row) + '   (3 jobs: measured | idle)')
    k, t, row = pairs[-1]
    report('at equal time (24 meas. ~ 40 us), idle minus measured, job 1', row[0][1] - row[0][0],
           'measuring costs more than idling at equal time, with the stated length')


# ================================================================ PART A: coherence cost per measurement, ibm_kingston

def king_class_curves(run):
    ks = KING[F_LADDER]; wt = KING['laid_waits_us']; out = {}
    for cls in (1, 2, 3):
        Q = [c for c in run['spectators'] if c['links_to_nearest_measured'] == cls]
        out[cls] = dict(n=len(Q),
                        meas=[median([c[F_SPEC_K][i] for c in Q]) for i in range(len(ks))],
                        laid=[median([c['hand_after_laid_wait'][i] for c in Q]) for i in range(len(wt))])
    return out


def laid_stretched(laid, wins=(10, 20, 30)):
    """fade curve h(t) = exp(-(t/tau)^b) fitted to the laid-time points at the stated windows (least squares in h)."""
    wt = KING['laid_waits_us']
    x = [t for t in wt if t in wins]; y = [laid[wt.index(t)] for t in x]
    best = None
    for ib in range(0, 271):
        b = 0.3 + 2.7 * ib / 270
        for it in range(0, 301):
            tau = math.exp(math.log(3) + (math.log(3000) - math.log(3)) * it / 300)
            c = sum((yi - math.exp(-(xi / tau) ** b)) ** 2 for xi, yi in zip(x, y))
            if best is None or c < best[0]: best = (c, tau, b)
    c, tau, b = best
    step = 0.1
    while step > 1e-4:
        moved = False
        for dt, db in ((1 + step, 1), (1 / (1 + step), 1), (1, 1 + step), (1, 1 / (1 + step))):
            t2, b2 = tau * dt, b * db
            if not (0.3 <= b2 <= 3.0 and 1 <= t2 <= 1e4): continue
            c2 = sum((yi - math.exp(-(xi / t2) ** b2)) ** 2 for xi, yi in zip(x, y))
            if c2 < c: c, tau, b, moved = c2, t2, b2, True
        if not moved: step /= 2
    return lambda t: math.exp(-(t / tau) ** b)


def score_run(run, k, ro_us, laid_model):
    """returns (kept ratio, lost ratio) averaged over the three distance classes."""
    ks = KING[F_LADDER]; wt = KING['laid_waits_us']; kr, lr = [], []
    for cur in king_class_curves(run).values():
        hk = loginterp(k, ks, cur['meas'])                       # k = 27 lies between the 24 and 32 points
        if laid_model == 'fit10-30':
            hl = laid_stretched(cur['laid'])(k * ro_us)
        elif laid_model == 'fit10-120':
            hl = laid_stretched(cur['laid'], (10, 20, 30, 60, 120))(k * ro_us)
        else:                                                    # 'measured': straight between the measured windows (1 at 0 us)
            hl = loginterp(k * ro_us, [0.0] + wt, [1.0] + cur['laid'])
        kr.append(hk / hl); lr.append((1 - hk) / (1 - hl))
    return mean(kr), mean(lr)


def part_a():
    print('\nA.  MEASUREMENTS vs THE SAME SPAN OF IDLE, ibm_kingston, six jobs (two days, three jobs each)')
    runs = KING['runs']
    print('  jobs: ' + '; '.join(f"{r['day']} run {r['run']} = {r['job_id']}" for r in runs))
    print('  each job: untouched qubits split by wiring distance to the nearest measured qubit (1, 2, 3+ links); per group the median\n'
          '  fraction of coherence kept after k measurements and after laid idle of 10/20/30/60/120 us; three groups averaged.')
    ro = {'2026-10-04': KING[F_MEAS_US]['2026-10-04'], '2026-10-06': KING[F_MEAS_US]['2026-10-06']}
    stated = KING[F_STATED_US]
    print(f"  length of one measurement: stated {stated} us; measured {ro['2026-10-04']} (4 Oct) and {ro['2026-10-06']} (6 Oct) (stored values, not recomputed here).")
    print('  There is NO 27-measurement point in these jobs (ladder 1,2,3,4,6,8,12,16,24,32); 27 is read between the 24 and 32 points.')
    table = {}
    for model, nm in (('fit10-30', 'idle curve fitted to the 10/20/30 us windows'),
                      ('fit10-120', 'idle curve fitted to all windows 10..120 us'),
                      ('measured', 'idle read straight off the measured windows (no fit)')):
        for lk in ('measured', 'stated'):
            rows = {}
            for k in (3, 9, 27):
                kv, lv = [], []
                for r in runs:
                    L = ro[r['day']] if lk == 'measured' else stated
                    a, b = score_run(r, k, L, model); kv.append(a); lv.append(b)
                rows[k] = (kv, lv)
            table[(model, lk)] = rows
    for (model, lk), rows in table.items():
        print(f"\n  idle reference: {model:<9} | measurement length: {lk}")
        for k in (3, 9, 27):
            kv, lv = rows[k]
            print(f"    k={k:>2}   kept(measured)/kept(idle): " + ' '.join(f'{v:.3f}' for v in kv) + f"   range {min(kv):.3f}-{max(kv):.3f}"
                  f"  |  lost/lost: " + ' '.join(f'{v:.3f}' for v in lv) + f"   range {min(lv):.3f}-{max(lv):.3f}")
    print('\n  Summary ratios (cost = measurement / same span of idle; idle fitted to 10-30 us, measured length):')
    base = table[('fit10-30', 'measured')]
    lo = min(min(base[3][0]), min(base[9][0])); hi = max(max(base[3][0]), max(base[9][0]))
    report('k <= 9, kept/kept, low end', lo)
    report('k <= 9, kept/kept, high end', hi)
    lostlo = min(min(base[3][1]), min(base[9][1])); losthi = max(max(base[3][1]), max(base[9][1]))
    report('k <= 9, lost/lost, low end', lostlo, 'k=3 alone is far above 1: loss is tiny at k=3, ratio unstable')
    report('k <= 9, lost/lost, high end', losthi)
    l9 = base[9][1]
    report('k = 9 only, lost/lost range, low end', min(l9))
    report('k = 9 only, lost/lost range, high end', max(l9))
    k27_kept = base[27][0]; k27_lost = base[27][1]
    report('k = 27, lost/lost, low end', min(k27_lost), 'six single jobs')
    report('k = 27, lost/lost, high end', max(k27_lost))
    d1 = [v for v, r in zip(k27_lost, runs) if r['day'] == '2026-10-04']; d2 = [v for v, r in zip(k27_lost, runs) if r['day'] == '2026-10-06']
    report('k = 27 lost/lost, mean of day 1 (3 jobs)', mean(d1))
    report('k = 27 lost/lost, mean of day 2 (3 jobs)', mean(d2))
    d1k = [v for v, r in zip(k27_kept, runs) if r['day'] == '2026-10-04']; d2k = [v for v, r in zip(k27_kept, runs) if r['day'] == '2026-10-06']
    report('k = 27 kept/kept, day 1 (3 jobs)', mean(d1k))
    report('k = 27 kept/kept, day 2 at its own measured length (3 jobs)', mean(d2k))
    report('k = 27, 1 / (kept ratio), mean of day 1', mean([1 / v for v, r in zip(k27_kept, runs) if r['day'] == '2026-10-04']))
    report('k = 27, 1 / (kept ratio), mean of day 2', mean([1 / v for v, r in zip(k27_kept, runs) if r['day'] == '2026-10-06']))
    free = table[('measured', 'measured')]
    report('k = 27, NO fitted curve (idle read off measured windows), kept/kept, low', min(free[27][0]), 'model-free reading')
    report('k = 27, NO fitted curve, kept/kept, high', max(free[27][0]))
    report('k = 27, NO fitted curve, lost/lost, mean of six', mean(free[27][1]))

    # direct, model-free equal-time check using only measured points: 60 us of idle vs the k that makes k*L = 60 us
    print('\n  Direct check with no fit: coherence kept after 60 us of laid idle vs after k = 60/L measurements (k read between measured points):')
    ks = KING[F_LADDER]
    for r in runs:
        L = ro[r['day']]; kk = 60.0 / L; row = []
        for cur in king_class_curves(r).values():
            row.append(loginterp(kk, ks, cur['meas']) / cur['laid'][KING['laid_waits_us'].index(60)])
        print(f"      {r['day']} run {r['run']}  k={kk:.1f}: measured/idle (1,2,3+ links) " + ' '.join(f'{v:.3f}' for v in row) + f'   mean {mean(row):.3f}')


def main():
    print('verify_cost.py  --  every figure below is recomputed from data/*.json')
    part_a()
    part_b()


if __name__ == '__main__':
    main()
