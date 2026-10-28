#!/usr/bin/env python3
"""verify_processors_apparatus.py -- recompute the processor-sort, decay-time, idle-control, non-exponential-decay,
gate-error and simulated-clock numbers from the compact JSON in data/.

    python3 verify_processors_apparatus.py

Standard library only. No account, no network, no numpy. Reads data/*.json (about 0.7 MB) and prints each
recomputed quantity with the processor, job and sample size it comes from, followed by short notes on what the
data do and do not show.
"""
import json, math, os, random, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
D = lambda f: json.load(open(os.path.join(HERE, 'data', f)))
READOUTS = [1, 2, 3, 4, 6, 9]


def note(s):
    print('   NOTE: ' + s)


def head(s):
    print('\n' + '=' * 100 + '\n' + s + '\n' + '=' * 100)


# ------------------------------------------------------------------------------------------------ small numerics
def linfit(x, y):
    n = len(x); mx = sum(x) / n; my = sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx


def nanmedian(v):
    v = [a for a in v if a == a]
    return st.median(v) if v else float('nan')


def unwrap(seq):
    out = [seq[0]]; corr = 0.0
    for a, b in zip(seq, seq[1:]):
        d = b - a; dd = (d + math.pi) % (2 * math.pi) - math.pi
        if dd == -math.pi and d > 0: dd = math.pi
        if abs(d) >= math.pi: corr += dd - d
        out.append(b + corr)
    return out


def quad_no_intercept(t, y):
    """least squares y = c0 t + c1 t^2 (quadratic through the origin, fitted to the free-evolution angles)."""
    s2 = sum(a ** 2 for a in t); s3 = sum(a ** 3 for a in t); s4 = sum(a ** 4 for a in t)
    b0 = sum(a * b for a, b in zip(t, y)); b1 = sum(a * a * b for a, b in zip(t, y))
    det = s2 * s4 - s3 * s3
    return (b0 * s4 - b1 * s3) / det, (s2 * b1 - s3 * b0) / det


def equiv_time(c0, c1, a):
    """real root of c1 x^2 + c0 x = a, inside (-5, 60) us, nearest to the linear estimate a/c0."""
    if c0 == 0: return float('nan')
    if c1 == 0: r = [a / c0]
    else:
        disc = c0 * c0 + 4 * c1 * a
        r = [] if disc < 0 else [(-c0 + s * math.sqrt(disc)) / (2 * c1) for s in (1, -1)]
    r = [x for x in r if -5 < x < 60]
    return min(r, key=lambda x: abs(x - a / c0)) if r else float('nan')


def spectator_readout_length(angles):
    """angles: [T10,T20,T30,K1,K2,K3,K4,K6,K9] radians -> readout length (us) by the slope over readouts 1..9; nan if <4 usable."""
    uw = unwrap([0.0] + angles[:3])
    c0, c1 = quad_no_intercept([10.0, 20.0, 30.0], uw[1:4])
    teq = [equiv_time(c0, c1, a) for a in angles[3:]]
    ok = [i for i, v in enumerate(teq) if v == v]
    if len(ok) < 4: return float('nan')
    return linfit([READOUTS[i] for i in ok], [teq[i] for i in ok])


def measured_qubit_readout_length(row):
    """row: [AW10,AW20,AW40, A1..A9 levels] of a measured qubit -> its own-clock readout length (us), or None."""
    clip = lambda v, lo, hi: min(max(v, lo), hi)
    L = [1.0] + [clip(v, 0.02, 1.0) for v in row[:3]]
    sl = linfit([0.0, 10.0, 20.0, 40.0], [math.log(v) for v in L])
    if sl >= 0: return None
    teq = [math.log(clip(v, 0.02, 1.0)) / sl for v in row[3:]]
    return linfit([float(k) for k in READOUTS], teq)


def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den


def ranks(v):                       # ordinal ranks, like argsort(argsort(.)) (no tie handling)
    order = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
    for k, i in enumerate(order): r[i] = k
    return r


def spearman(a, b):
    """with average ranks for ties"""
    def avg(v):
        s = sorted(v); out = []
        for x in v:
            idx = [i for i, y in enumerate(s) if y == x]; out.append(sum(idx) / len(idx))
        return out
    return pearson(avg(a), avg(b))


# ================================================================================================ 1. PROCESSORS SORT
# Median T2:T1 per chip from the chip calibration record (the daily calibration data of the processor vendor;
# the chip records themselves are not re-derived here).
T2T1 = {'ibm_marrakesh': 0.519, 'ibm_kingston': 0.638, 'ibm_fez': 0.683, 'ibm_phoenix': 0.851, 'ibm_pittsburgh': 1.160, 'ibm_boston': 1.183}
# The 8 chip-days (5 chips) over which the three spectator-selection rules are compared by rank agreement.
RANK_SET = sorted([('ibm_boston', '4 Oct'), ('ibm_fez', '4 Oct'), ('ibm_fez', '6 Oct'), ('ibm_kingston', '4 Oct'), ('ibm_kingston', '6 Oct'),
                   ('ibm_marrakesh', '4 Oct'), ('ibm_pittsburgh', '4 Oct'), ('ibm_pittsburgh', '5 Oct')])


def sort_section():
    head('1. PROCESSORS SORT   (ratio = measured qubit\'s own-clock readout length over the unmeasured spectators\' readout length)')
    data = D('spectator_and_measured_qubit_sort.json')
    rows = {}
    for cd in data['chipdays']:
        Ts, Wa, Wb, Ws = [], [], [], []
        for run in cd['runs']:
            spec_key = next(k for k in run if k.endswith('_rows') and k != 'far_rows')   # spectator rows
            W = [(r[1], spectator_readout_length(r[2:])) for r in run[spec_key]]
            Wa.append(nanmedian([w for _, w in W]))
            b = [w for n, w in W if n == 0]; s = [w for n, w in W if n >= 2]
            Wb.append(nanmedian(b) if b else float('nan')); Ws.append(nanmedian(s) if s else float('nan'))
            tt = [measured_qubit_readout_length(r[1:]) for r in run['far_rows']]; Ts.append(st.median([t for t in tt if t is not None]))
        T = st.mean(Ts)
        rows[(cd['chip'], cd['day'])] = (T / st.mean(Wa), T / st.mean(Wb) if all(x == x for x in Wb) else None,
                                         T / st.mean(Ws) if all(x == x for x in Ws) else None, [r['run'] for r in cd['runs']])
    print('   chip-day              runs          ratio, all spectators   ratio, no measured neighbour   ratio, 2+ measured neighbours')
    for k in sorted(rows):
        p, b, s, runs = rows[k]
        print('   %-14s %-6s %-14s %8.3f   %8s   %9s' % (k[0], k[1], ','.join(map(str, runs)), p, '%.3f' % b if b else 'n/a (layout not Heron)',
              '%.3f' % s if s else 'n/a'))
    print()
    # rank agreement across the 8 chip-days: ordinal ranks, Pearson
    P = [rows[k][0] for k in RANK_SET]; B = [rows[k][1] for k in RANK_SET]; S = [rows[k][2] for k in RANK_SET]
    print('   rank agreement, all spectators vs no measured neighbour (8 chip-days, 5 chips): %+.2f' % pearson(ranks(P), ranks(B)))
    print('   rank agreement, all spectators vs 2+ measured neighbours (8 chip-days, 5 chips): %+.2f' % pearson(ranks(P), ranks(S)))
    # per chip (mean of days)
    chip = {}
    for (c, d), v in rows.items(): chip.setdefault(c, []).append(v[0])
    print('\n   per-chip ratio, all spectators (every chip-day listed above):')
    for c in sorted(chip, key=lambda c: st.mean(chip[c])):
        print('     %-16s %s   mean %.3f   %s' % (c, ' '.join('%.3f' % x for x in chip[c]), st.mean(chip[c]),
              'calibration T2:T1 %.3f' % T2T1[c] if c in T2T1 else 'calibration T2:T1 not available'))
    lowg = ['ibm_boston', 'ibm_pittsburgh', 'ibm_phoenix']; highg = ['ibm_kingston', 'ibm_fez', 'ibm_marrakesh']
    lo = [x for c in lowg for x in chip[c]]; hi = [x for c in highg for x in chip[c]]
    print('   lowest ratio of Boston, Pittsburgh, Phoenix: %.3f' % min(lo))
    print('   highest ratio of Boston, Pittsburgh, Phoenix: %.3f' % max(lo))
    print('   lowest ratio of Kingston, Fez, Marrakesh: %.3f' % min(hi))
    print('   highest ratio of Kingston, Fez, Marrakesh: %.3f' % max(hi))
    print()
    print('   two groups or a continuum?')
    allv = sorted((x, c) for c in chip for x in chip[c])
    gaps = [(allv[i + 1][0] - allv[i][0], allv[i][1], allv[i + 1][1]) for i in range(len(allv) - 1)]
    print('     all chip-day values, sorted: ' + ' '.join('%.3f' % v for v, _ in allv))
    big = max(gaps)
    bnd = min(hi) - max(lo)
    print('     gap between the groups (highest of Boston/Pittsburgh/Phoenix to lowest of Kingston/Fez/Marrakesh): %.3f' % bnd)
    print('     largest gap anywhere: %.3f (between %s and %s)' % (big[0], big[1][4:], big[2][4:]))
    fez = chip['ibm_fez']
    print('     Fez alone moves %.3f between its two days, %.1f times the gap between the groups' % (max(fez) - min(fez), (max(fez) - min(fez)) / bnd))
    print('     spread inside the low group %.3f, inside the high group %.3f (with Aachen in the low group: %.3f)' % (max(lo) - min(lo), max(hi) - min(hi), max(lo) - min(chip['ibm_aachen'] + lo)))
    m6 = [c for c in T2T1]
    xs = [st.mean(chip[c]) for c in m6]; ys = [T2T1[c] for c in m6]
    rho = spearman(xs, ys)
    print('     rank correlation, chip mean ratio against calibration T2:T1, six chips: %+.3f' % rho)
    pred = sorted(m6, key=lambda c: -T2T1[c]); obs = sorted(m6, key=lambda c: st.mean(chip[c]))
    print('       order by T2:T1 high->low: ' + ' '.join(c[4:] for c in pred))
    print('       order by ratio low->high: ' + ' '.join(c[4:] for c in obs))
    # chance of the two-group split: label chips with T2:T1>1 (2 chips) vs the rest; how often would a random label put both in the lowest two?
    n = len(m6); tot = math.comb(n, 2)
    print('     chance that the 2 chips with T2:T1 > 1 are the 2 lowest in ratio by luck, among 6 chips: 1 in %d' % tot)
    ph = chip['ibm_phoenix']
    note('Aachen is the 7th chip and is computed above; its calibration T2:T1 is not available, so it cannot be placed on the T2:T1 axis.')
    note('Phoenix has T2:T1 = %.3f (below 1, like Kingston/Fez/Marrakesh) but its ratio (%.2f-%.2f) sits in the low group with Boston and Pittsburgh.'
         % (T2T1['ibm_phoenix'], min(ph), max(ph)))
    note('Within the high group the order by ratio is Marrakesh < Kingston < Fez, the reverse of the T2:T1 order Fez > Kingston > Marrakesh. The data support a split at T2:T1 = 1 on heavy-hex chips, not an ordering.')
    note('The two rank agreements above are between three ways of choosing spectators, over 8 chip-days of 5 chips; they are not agreement between the ratio and T2:T1. Phoenix and Aachen are not in them.')


# ================================================================================================ 2 & 4. TIME RUNS
def tau_of(l300):
    return -300.0 / math.log(min(max(l300, 0.02), 0.98))


def time_section():
    head('2. T1 DRIFT   (reference decay: fading time of a qubit left in |1> and read after 300 us, three copies a/b/c per job)')
    data = D('time_runs_ruler_copies.json')['runs']
    by = {}
    for run in data:
        taus = []
        for ci in range(3):
            base = 1 + 4 * ci         # row = [q, (30,100,300,600) x copy a, b, c, idle300, idle600]
            taus.append(st.median([tau_of(r[base + 2]) for r in run['rows']]))
        by.setdefault(run['chip'], []).append((run['run'], run['start'][:16], taus))
    for chip in ('ibm_marrakesh', 'ibm_kingston', 'ibm_pittsburgh'):
        runs = by[chip]
        print('\n   %s   fading time (us), median over readable qubits' % chip)
        print('     %-6s %-17s %8s %8s %8s   | mean of copies   range of copies (inside one job)' % ('run', 'started', 'copy a', 'copy b', 'copy c'))
        for r, t0, t in runs:
            print('     %-6d %-17s %8.1f %8.1f %8.1f   | %8.1f         %5.1f' % (r, t0, t[0], t[1], t[2], st.mean(t), max(t) - min(t)))
        copy_means = [st.mean([r[2][c] for r in runs]) for c in range(3)]
        read_means = [st.mean(r[2]) for r in runs]
        print('     mean over the %d jobs, copy a/b/c: %s      |  mean over copies, job by job: %s' % (len(runs), ' / '.join('%.1f' % x for x in copy_means), ' / '.join('%.1f' % x for x in read_means)))
        by[chip] = (runs, copy_means, read_means)
    print()
    m = by['ibm_marrakesh']; k = by['ibm_kingston']; p = by['ibm_pittsburgh']
    for i, c in enumerate(m[1]):
        print('   marrakesh: mean of copy %s over 3 jobs: %.1f us' % ('abc'[i], c))
    for i, c in enumerate(m[2]):
        print('   marrakesh job %d, mean over its three copies: %.1f us' % (i + 1, c))
    print('   kingston lowest copy mean: %.1f us' % min(k[1]))
    print('   kingston highest copy mean: %.1f us' % max(k[1]))
    pr = p[0]
    first = pr[0][2]
    print('   pittsburgh run %d, lowest of the three copies: %.1f us' % (pr[0][0], min(first)))
    print('   pittsburgh run %d, highest of the three copies: %.1f us' % (pr[0][0], max(first)))
    spread = lambda runs, c: max(r[2][c] for r in runs) - min(r[2][c] for r in runs)
    print('\n   how far each copy moves BETWEEN jobs (max - min of that copy over the jobs):')
    for nm, bb in (('marrakesh', m[0]), ('kingston', k[0]), ('pittsburgh', p[0])):
        print('     %-11s copy a/b/c: %s us' % (nm, ' / '.join('%.1f' % spread(bb, c) for c in range(3))))
    print('   largest job-to-job move of any copy, 3 chips: %.1f us' % max(spread(bb, c) for bb in (m[0], k[0], p[0]) for c in range(3)))
    within = [max(r[2]) - min(r[2]) for bb in (m[0], k[0], p[0]) for r in bb]
    print('   range of the three copies INSIDE a single job: min %.1f, max %.1f us (marrakesh jobs: %s)' % (min(within), max(within), ' / '.join('%.1f' % (max(r[2]) - min(r[2])) for r in m[0])))
    note('The three copies (a, b, c) are measured inside each job. Their means are averaged over jobs for Marrakesh and Kingston; for Pittsburgh the first job alone is quoted above.')
    note('The copy-to-copy difference inside one job is as large as the job-to-job movement (Marrakesh job 1: %.1f us).' % (max(m[0][0][2]) - min(m[0][0][2])))
    note('Pittsburgh has three jobs (runs 361-363); their copy-means are %s us.' % ' / '.join('%.1f' % x for x in p[2]))
    return data


def bend_section(data):
    head('4. NON-EXPONENTIAL DECAY   (|1> level at 30 us below that qubit-copy\'s own 300 us single-exponential extrapolation)')
    first_set = {('ibm_marrakesh', 301), ('ibm_marrakesh', 302), ('ibm_marrakesh', 303), ('ibm_kingston', 301), ('ibm_kingston', 302), ('ibm_kingston', 303),
                 ('ibm_pittsburgh', 361)}
    shares = {}
    print('   chip           run   readable qubit-copies   share below at 30 us')
    for run in data:
        n = b = 0
        for r in run['rows']:
            for ci in range(3):
                base = 1 + 4 * ci; l30, l300 = r[base], r[base + 2]
                tau = tau_of(l300); pred = math.exp(-30.0 / tau)
                n += 1; b += (l30 - pred) < 0
        shares[(run['chip'], run['run'])] = 100.0 * b / n
        print('   %-14s %-5d %-22d %5.1f %%' % (run['chip'], run['run'], n, shares[(run['chip'], run['run'])]))
    print()
    six = [shares[k] for k in sorted(first_set)]
    print('   marrakesh runs 301-303, kingston runs 301-303, pittsburgh run 361 (7 runs): lowest %.0f %%, highest %.0f %%' % (min(six), max(six)))
    others = {k: v for k, v in shares.items() if k not in first_set}
    print('   the other time-design runs on disk: lowest %.0f %%, highest %.0f %% over %d runs' % (min(others.values()), max(others.values()), len(others)))
    pit = [v for k, v in shares.items() if k[0] == 'ibm_pittsburgh']
    print('   Pittsburgh, all three jobs: %s %%' % ' / '.join('%.0f' % x for x in pit))
    note('The share below the extrapolation never drops under 50 % in any run, so the sign of the effect is stable; the exact range depends on which jobs are counted.')
    note('Idle |0> after 300 us (the "left alone" arm of the same jobs) is in section 3.')


# ================================================================================================ 3. PARK + IDLE
def book(p300, nseg):
    q = p300 ** (1.0 / nseg); on, off = 1.0, 0.0
    for _ in range(nseg): on, off = off + on * (1 - q), on * q
    return on


def idle_section(pdata):
    head('3. IDLE CONTROL')
    print('   (a) alternating-state runs 71-73, ibm_pittsburgh: qubit starts in |1>, an X gate after each interval so it spends alternate intervals in |0>,')
    print('       300 us in all; measured level minus that qubit\'s own decay bookkeeping (its own 300-us nothing-done level, cut into the intervals).')
    arms = {'P1': ('a  parked once       (2 x 150 us)', 2), 'P3': ('b  alternate quarters (4 x 75 us)', 4), 'P9': ('c  alternate tenths   (10 x 30 us)', 10)}
    for key, (label, n) in arms.items():
        means, widths = [], []
        for run in pdata['park']:
            sep = run['sep']; g = [i for i, s in enumerate(sep) if s > 4]
            dl, sq = [], []
            for i in g:
                p300 = min(max(run['lv_N0'][i], 0.02), 0.98)
                d = run['lv_' + key][i] - book(p300, n); dl.append(d)
                eps = 1e-4; dd = (book(min(p300 + eps, 0.98), n) - book(max(p300 - eps, 0.02), n)) / (2 * eps)
                sq.append(run['se_' + key][i] ** 2 + (dd * run['se_N0'][i]) ** 2)
            means.append(st.mean(dl)); widths.append(math.sqrt(sum(sq)) / len(g))
        m = st.mean(means); w = math.sqrt(sum(x * x for x in widths)) / 3
        print('     %s: reads %s (readable %s)' % (label, ' '.join('%+.4f' % x for x in means), '/'.join(str(sum(1 for s in r['sep'] if s > 4)) for r in pdata['park'])))
        print('        residual, three-read mean: %+.4f' % m)
        print('        width (analytic shot error): %.4f' % w)
    note('These runs prepare |1> and flip it with X gates, so the |0> time is every other interval; no run here prepares |0> and leaves it. One chip (Pittsburgh), one day, 3 jobs.')
    print()
    print('   (b) idle |0> for 300 us: the "left alone 300 us" arm of run 121, ibm_pittsburgh (one job).')
    rows = pdata['idle_run121']['rows']
    z = [r[1] for r in rows]; se = [math.hypot(r[2], r[3]) for r in rows]   # level error and reference error
    mean = st.mean(z); w = math.sqrt(sum(x * x for x in se)) / len(z)
    print('     idle |0> after 300 us, gain in percent: %+.3f' % (100 * mean))
    print('     width in percent: %.3f' % (100 * w))
    print('     readable qubits: %d' % len(z))
    # the same arm in the time-design jobs on the other chips
    td = D('time_runs_ruler_copies.json')['runs']
    print('   the same idle-|0> arm (300 us) in the time-design jobs on all chips, mean over readable qubits (percent):')
    per = {}
    for run in td:
        per.setdefault(run['chip'], []).append(100 * st.mean([r[13] for r in run['rows']]))
    for c in sorted(per): print('     %-15s %s' % (c, ' '.join('%+.2f' % x for x in per[c])))
    allj = [x for v in per.values() for x in v]
    note('Run 121 is one Pittsburgh job. In the %d time-design jobs above the same arm reads %+.2f to %+.2f %% (positive in %d).' % (len(allj), min(allj), max(allj), sum(1 for x in allj if x > 0)))


# ================================================================================================ 5. GATE ERROR
def gate_section(pdata):
    head('5. GATE ERROR GROWING WITH GATE COUNT   (push run 91-93, ibm_pittsburgh, back-to-back X gates, no wait)')
    acc = {}
    for run in pdata['push']:
        rows = run['rows']; names = run['columns']
        idx = {n: i for i, n in enumerate(names)}
        for arm, n in (('OFF', 2), ('OFF', 10), ('ON', 2), ('ON', 10)):
            a = 'lv %s|p%d' % (arm, n); b = 'lv %s|p0' % arm
            sa = 'se %s|p%d' % (arm, n); sb = 'se %s|p0' % arm
            d = [r[idx[a]] - r[idx[b]] for r in rows]
            w = math.sqrt(sum(r[idx[sa]] ** 2 + r[idx[sb]] ** 2 for r in rows)) / len(rows)
            acc.setdefault((arm, n), []).append((st.mean(d), w))
    val = {}
    for key, lst in acc.items():
        m = st.mean(x for x, _ in lst); w = math.sqrt(sum(x * x for _, x in lst)) / len(lst)
        val[key] = (m, w)
    for key in (('OFF', 2), ('ON', 2), ('OFF', 10), ('ON', 10)):
        print('   shift after %2d back-to-back X gates, start %s (percent): %+.2f' % (key[1], key[0], 100 * val[key][0]))
    # weighted least squares of o_N = N^2 e2/4 + N b  (e2 = eps^2) on the four numbers
    pts = [(2, val[('OFF', 2)]), (2, (-val[('ON', 2)][0], val[('ON', 2)][1])), (10, val[('OFF', 10)]), (10, (-val[('ON', 10)][0], val[('ON', 10)][1]))]
    A11 = A12 = A22 = B1 = B2 = 0.0
    for N, (y, s) in pts:
        w = 1 / s ** 2; f1 = N * N / 4.0; f2 = float(N)
        A11 += w * f1 * f1; A12 += w * f1 * f2; A22 += w * f2 * f2; B1 += w * f1 * y; B2 += w * f2 * y
    det = A11 * A22 - A12 * A12
    e2 = (B1 * A22 - B2 * A12) / det; b = (A11 * B2 - A12 * B1) / det
    se_e2 = math.sqrt(A22 / det); se_b = math.sqrt(A11 / det)
    eps = math.sqrt(e2); eps_deg = math.degrees(eps); se_deg = math.degrees(se_e2 / (2 * eps))
    print('   rotation shortfall per gate (degrees): %.2f' % eps_deg)
    print('     (+-%.2f deg)' % se_deg)
    print('     count-linear term b in  shift(N) = N^2 eps^2/4 + N b :  b = %.3f +- %.3f %%' % (100 * b, 100 * se_b))
    print('   b x 1 (as population error per gate), percent: %.3f' % (100 * b))
    print('   4 x b (percent): %.3f' % (400 * b))
    print('      its width (4 x width of b), percent: %.3f' % (400 * se_b))
    print('     upper edge at two widths, 4 x (b + 2 se): %.2f %%' % (400 * (b + 2 * se_b)))
    note('4 x b is a measured value, not an upper bound; the upper edge at two widths is given above. It covers EVERY count-growing error per gate (leakage, fading during the pulse, pulse noise), not leakage alone.')
    note('The fit assumes the |2> population after N gates is the same for an OFF and an ON start; one chip (Pittsburgh).')


# ================================================================================================ 6. SIMULATED CLOCKS
def simulate(plant, noise, seed, n_clocks=100, repeats=30, O=2.18):
    rng = random.Random(seed)
    fit_all, w39_all = [], []
    for rep in range(repeats):
        wfit, w39 = [], []
        for _ in range(n_clocks):
            rate = math.radians(rng.uniform(2.3, 5.6))               # rad per us
            ang = lambda t: rate * t + rng.gauss(0, noise)
            T = [ang(t) for t in (10.0, 20.0, 30.0)]
            K = [ang(k * O * (1.0 + plant)) for k in READOUTS]
            a = unwrap([0.0] + T)
            c0, c1 = quad_no_intercept([10.0, 20.0, 30.0], a[1:4])
            teq = [equiv_time(c0, c1, x) for x in K]
            ok = [i for i, v in enumerate(teq) if v == v]
            if len(ok) >= 4: wfit.append(linfit([READOUTS[i] for i in ok], [teq[i] for i in ok]))
            if teq[2] == teq[2] and teq[5] == teq[5]: w39.append((teq[5] - teq[2]) / 6.0)
        fit_all.append(nanmedian(wfit) / O); w39_all.append(nanmedian(w39) / O)
    return st.mean(fit_all), st.pstdev(fit_all), st.mean(w39_all), st.pstdev(w39_all)


def simulation_section():
    head('6. SIMULATED CLOCKS   (synthetic spectator clocks read with the same inversion as the real data)')
    print('   100 synthetic spectator clocks, rates uniform 2.3-5.6 deg/us, readouts of exactly 2.18 us x (1 + planted), 30 repeats per seed, 8 seeds,')
    print('   read with the same inversion as the real data (quadratic through the origin on 10/20/30 us, roots per readout, slope over 1..9 readouts).')
    print('   Gaussian phase noise added to every measured angle, at 0.06 rad and at 0.12 rad.')
    print('   %-38s %-34s %-34s' % ('case', 'fit 1-9: mean over seeds [min,max]', '3-and-9: mean over seeds [min,max]'))
    for label, plant, noise in (('nothing planted, no noise', 0.0, 0.0), ('nothing planted, phase noise 0.06 rad', 0.0, 0.06), ('nothing planted, phase noise 0.12 rad', 0.0, 0.12),
                                ('4 % planted, no noise', 0.04, 0.0), ('4 % planted, phase noise 0.06 rad', 0.04, 0.06), ('4 % planted, phase noise 0.12 rad', 0.04, 0.12)):
        rs = [simulate(plant, noise, 100 + s) for s in range(8)]
        f = [r[0] for r in rs]; g = [r[2] for r in rs]
        print('   %-38s %.4f  [%.4f, %.4f]          %.4f  [%.4f, %.4f]' % (label, st.mean(f), min(f), max(f), st.mean(g), min(g), max(g)))
    note('With NO noise a planted 4 % reads back 1.0400 exactly; with noise it reads a few tenths of a percent low, and the unplanted case a few tenths low of 1.000. The seed-to-seed range is about +-0.003.')
    note('This tests the spectator-clock readout-length estimator only. It does not test the ratio of section 1.')


# ================================================================================================
def main():
    print(__doc__.split('\n\n')[0])
    sort_section()
    td = time_section()
    pd = D('pittsburgh_park_idle_push.json')
    idle_section(pd)
    bend_section(td)
    gate_section(pd)
    simulation_section()


if __name__ == '__main__':
    main()
