#!/usr/bin/env python3
"""verify_gate_cost.py -- recompute the gate-cost numbers from the released data.

    python3 verify_gate_cost.py

No account, no network, no dependencies beyond the standard library. Reads data/gate_cost.json
(per-qubit levels from 17 hardware jobs, ~300 KB) and recomputes each figure.

Estimator (plain words). A qubit starts in its ground state (population 0). We apply n X gates
(n = 2, 6, 10), waiting 30 microseconds after each. Each X gate flips ground<->excited, and during
every 30 us wait the excited population decays by a factor q (the survival over 30 us). Predicted
excited population after n gates = recursion  on <- off + on*(1-q), off <- on*q  (start on=0, off=1),
i.e. the gates are assumed perfect and only decay acts. The RESIDUAL is measured minus predicted,
averaged over qubits. The SLOPE PER GATE is the least-squares slope of the residual against n.
q comes from a reference measurement in the same job (qubit prepared excited and left alone):
  "own interval"  : q = the 30 us survival measured directly
  "90 us ruler"   : q = (survival after 90 us) ** (30/90)
  "one exponential": q = (survival after 300 us) ** (30/300)
Survival values are clipped to [0.02, 1]; for ibm_pittsburgh on the own-interval ruler the upper clip is 1.2.
Widths are 1-sigma, propagated from shot noise on every level (no resampling; first order).
"""
import json, math, os, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
NS = (2, 6, 10)



def book(q, n):
    on, off = 0.0, 1.0
    for _ in range(n):
        on, off = off + on * (1 - q), on * q
    return on


def clip(x, lo, hi): return min(max(x, lo), hi)


def qfun(ruler, r, i):
    hi = 1.2 if r['chip'] == 'ibm_pittsburgh' and ruler == 'own' else 1.0
    if ruler == 'own':  return clip(r['R30'][i], 0.02, hi), r['se_R30'][i], 1.0
    if ruler == 'r90':  return clip(r['R90'][i], 0.02, 1.0) ** (1 / 3), r['se_R90'][i], 1 / 3
    if ruler == 'r300': return clip(r['R300'][i], 0.02, 1.0) ** (1 / 10), r['se_R300'][i], 1 / 10


def residuals(r, ruler):
    """per-qubit residual at n=2,6,10 and per-qubit variance of (slope contribution)."""
    rows = []
    for i in range(len(r['qubits'])):
        q, seq, p = qfun(ruler, r, i)
        # derivative of q wrt the raw survival s: p * s^(p-1) (where not clipped) -> numeric
        d = []
        for n in NS:
            d.append(r['N%d' % n][i] - book(q, n))
        h = 1e-5
        dq = [(book(min(q + h, 1.0), n) - book(q - h, n)) / (min(q + h, 1.0) - (q - h)) for n in NS]
        # chain rule: q = s**p  ->  dq/ds = p*s**(p-1)
        s_raw = {'own': r['R30'][i], 'r90': r['R90'][i], 'r300': r['R300'][i]}[ruler]
        dqds = p * clip(s_raw, 0.02, 1.0) ** (p - 1) if p != 1.0 else 1.0
        rows.append((d, dq, dqds, seq))
    return rows


def slope_stats(r, ruler):
    rows = residuals(r, ruler); N = len(rows)
    dmean = [sum(row[0][k] for row in rows) / N for k in range(3)]
    slope = (dmean[2] - dmean[0]) / 8.0            # least squares through n = 2, 6, 10
    var = 0.0
    for i, (d, dq, dqds, seq) in enumerate(rows):
        var += (r['se_N10'][i] ** 2 + r['se_N2'][i] ** 2 + ((dq[2] - dq[0]) * dqds * seq) ** 2) / 64.0
    width = math.sqrt(var) / N
    per_bit = []
    for d, *_ in rows:
        per_bit.append((d[2] - d[0]) / 8.0)
    return slope, width, dmean, per_bit


def trim_mean(x, cut=0.2):
    x = sorted(x); k = int(cut * len(x)); return st.mean(x[k:len(x) - k])


def avg_slopes(r, ruler='own'):
    rows = residuals(r, ruler); cols = [[row[0][k] for row in rows] for k in range(3)]
    f = lambda agg: (agg(cols[2]) - agg(cols[0])) / 8.0
    return f(st.mean), f(st.median), f(trim_mean)


def binom_two_sided(k, n):
    """exact two-sided sign test p-value against p = 0.5."""
    kk = min(k, n - k); tail = sum(math.comb(n, j) for j in range(kk + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def main():
    D = json.load(open(os.path.join(HERE, 'data', 'gate_cost.json')))['runs']
    by = {(r['chip'], r['run']): r for r in D}
    S = {}
    for r in D:
        for ruler in ('own', 'r90', 'r300'):
            S[(r['chip'], r['run'], ruler)] = slope_stats(r, ruler)
    g = lambda chip, run, ruler='own': S[(chip, run, ruler)]

    def line(label, got, wid=None, fmt='%+.5f'):
        s = '  %-44s %s' % (label, fmt % got)
        if wid is not None: s += ' (shot-noise width %.5f)' % wid
        print(s)

    print('GATE COST: slope of the residual excited population per X gate (own-interval estimator unless stated)')
    print('1. Slopes per gate. Widths are first-order shot noise on the levels only; they do not resample the reference readouts.')
    s, w, *_ = g('ibm_kingston', 392);      line('kingston day 1 (run 392)', s, wid=w)
    ks = [g('ibm_kingston', k) for k in (721, 722, 723)]
    m = st.mean(x[0] for x in ks); w = math.sqrt(sum(x[1] ** 2 for x in ks)) / 3
    line('kingston day 2, mean of runs 721-723', m, wid=w)
    s, w, *_ = g('ibm_fez', 371);           line('fez (run 371)', s, wid=w)
    for k, key in zip((711, 712, 713), ('boston read 1 (day 2)', 'boston read 2 (day 2)', 'boston read 3 (day 2)')):
        s, w, *_ = g('ibm_boston', k); line('%s (run %d)' % (key, k), s, wid=w)
    px = [g('ibm_phoenix', k) for k in (751, 752, 753)]
    pm_own = st.mean(x[0] for x in px); pm_90 = st.mean(g('ibm_phoenix', k, 'r90')[0] for k in (751, 752, 753))
    line('phoenix day 2, 3-read mean, own interval', pm_own)
    line('phoenix day 2, 3-read mean, 90 us ruler', pm_90)
    s, w, *_ = g('ibm_pittsburgh', 121);    line('pittsburgh 2 Oct (run 121)', s, wid=w)

    print('\n2. Further slopes')
    for chip, run in (('ibm_boston', 441), ('ibm_phoenix', 531)):
        s, w, *_ = g(chip, run)
        print('        %-13s day 1 run %d own interval  %+.5f +- %.5f  (%.1f sigma)' % (chip, run, s, w, abs(s) / w))
    for chip, run in (('ibm_phoenix', 531), ('ibm_phoenix', 751), ('ibm_phoenix', 752), ('ibm_phoenix', 753), ('ibm_kingston', 392), ('ibm_boston', 441)):
        print('        %-13s run %d  own %+.5f   90us %+.5f   one-exponential %+.5f' % ((chip, run) + tuple(g(chip, run, k)[0] for k in ('own', 'r90', 'r300'))))
    pm300 = st.mean(g('ibm_phoenix', k, 'r300')[0] for k in (751, 752, 753))
    print('        phoenix day-2 3-read mean on the one-exponential ruler: %+.5f' % pm300)

    print('\n3. Sign of the slope across estimators, reads and heavy-hex chips')
    core = [(c, k) for c, ks_ in (('ibm_fez', [371]), ('ibm_kingston', [392, 721, 722, 723]), ('ibm_boston', [441, 711, 712, 713])) for k in ks_]
    heavy = core + [('ibm_aachen', k) for k in (841, 842, 843)]
    def neg_report(runs, name):
        allneg = 0; ruler_neg = 0
        for c, k in runs:
            a = avg_slopes(by[(c, k)])
            rn = all(g(c, k, x)[0] < 0 for x in ('own', 'r90', 'r300')); an = all(x < 0 for x in a)
            ruler_neg += rn; allneg += rn and an
        print('  %s: %d of %d runs negative on all 3 rulers; %d of %d also negative on mean, median and 20%%-trim over qubits' % (name, ruler_neg, len(runs), allneg, len(runs)))
        for c, k in runs:
            a = avg_slopes(by[(c, k)]); bad = [x for x in ('own', 'r90', 'r300') if g(c, k, x)[0] >= 0] + [n for n, x in zip(('mean', 'median', 'trim'), a) if x >= 0]
            if bad: print('        not negative: %s %d on %s' % (c, k, bad))
    neg_report(core, 'fez, kingston, boston')
    neg_report(heavy, 'fez, kingston, boston and aachen (3 reads)')
    print('     single reads below 2 sigma on the own-interval ruler (shot-noise width):', ', '.join('%s %d (%.1f sigma)' % (c.replace('ibm_', ''), k, abs(g(c, k)[0]) / g(c, k)[1]) for c, k in heavy if abs(g(c, k)[0]) / g(c, k)[1] < 2) or 'none')

    print('\n4. Share of individual qubits with a negative own slope')
    rng_ = {}
    for c, k in heavy + [('ibm_phoenix', 531), ('ibm_phoenix', 751), ('ibm_phoenix', 752), ('ibm_phoenix', 753)]:
        pb = g(c, k)[3]; n = len(pb); neg = sum(1 for x in pb if x < 0)
        rng_[(c, k)] = (100.0 * neg / n, binom_two_sided(neg, n))
        print('        %-13s run %4d  %3d qubits  %3.0f%% negative  (sign-test p = %.3f)' % (c, k, n, rng_[(c, k)][0], rng_[(c, k)][1]))
    core_chips = [(c, k) for c, k in heavy if c != 'ibm_aachen']
    lo = min(rng_[x][0] for x in core_chips); hi = max(rng_[x][0] for x in core_chips)
    lo_all = min(rng_[x][0] for x in heavy); hi_all = max(rng_[x][0] for x in heavy)
    print('  range over fez, kingston, boston: %.0f-%.0f %%; over all heavy-hex runs incl. aachen: %.0f-%.0f %%' % (lo, hi, lo_all, hi_all))
    pl = [rng_[x][0] for x in rng_ if x[0] == 'ibm_phoenix']
    print('     phoenix runs fall at %.0f-%.0f %%' % (min(pl), max(pl)))

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
