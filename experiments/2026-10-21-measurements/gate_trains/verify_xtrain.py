#!/usr/bin/env python3
"""verify_xtrain.py -- recompute every number in main.tex Sec. "Spaced X-gate trains and the decay-bookkeeping control" (sec:xtrain),
the 2 October idle-level sentence of Sec. "Idle decay: shape ..." (paragraph Shape), the equator-start sentence of Sec. "Gate cost ..."
(paragraph Back-to-back gates) and the abstract's "over-correcting by 1 to 2% on marrakesh", from the raw IBM job results.

    /Users/heuremen/bones-forge/qvenv/bin/python3 verify_xtrain.py        (numpy only; ~1-2 minutes)

Writes results_xtrain.text next to this file and prints it. Reads only ./data/ (a verbatim copy of the raw per-circuit IQ results, run
meta and job records; see data/<chip>/run<N>_{results.npz,meta.json,job.json}).

RAW RECORD. results.npz holds one array per circuit: shots x qubits of complex IQ points (key 'OFF|p0__L1' = start |0>, no gate, one look).
LEVELS. For each qubit the OFF|p0 cloud is level 0 and the ON|p0 cloud level 1 (same job): level = Re[(x - c0) conj(d)] / |d|^2 with
c0 = mean of OFF|p0, d = mean(ON|p0) - c0. A qubit is READABLE when |d| exceeds 4 widths of the OFF|p0 cloud (width = sqrt(mean|x-c0|^2 / 2)).
All means are over readable qubits. (Run 121 uses its own frame circuits F0 / F1 with the same rule.)
BOOKKEEPING (the conventions the notebook / reader scripts state, parallax/out/qubit_fifth_run/{readclose,readtime,readnoname}.py):
  ONE   one exponential through the qubit's own 300 us start-|1> level (the 'ON left alone 300 us' circuit), q(spacing) = L300 ** (spacing/300)
        (L300 clipped to [0.02, 0.98]); the train is then the recursion  on <- off + on (1-q), off <- on q  over ten intervals, started from the
        qubit's start (0 for OFF, 1 for ON).
  CURVE the job's own start-|1> decay curve, read at 0, 30, 100, 300 and 600 us (one circuit each; level 1 at 0 us), log-linear between points,
        q(spacing) read off it (levels clipped to [0.02, 1.0]); below 30 us it is interpolated between the 0 and 30 us points.
  Three ruler copies (a, b, c) of the decay curve exist in the 4 October jobs; each bookkeeping is built from each copy separately, and the
  copies are averaged.
WIDTHS are 1 sigma from resampling the shots of every circuit (fixed seed, NB resamples); three jobs combine as sqrt(sum w^2)/3 (the reader
scripts' rule). The reference scripts used fewer resamples (30 to 100) and an unseeded-in-this-file generator, so a printed width is expected to
agree only to a few percent; a width is marked 'ok(boot)' when it agrees within 15%.
"""
import json, math, os, sys, warnings
import numpy as np

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
NB = 200
RNG = np.random.default_rng(20261022)
OUT = []


def P(s=''):
    OUT.append(s); print(s)


# ------------------------------------------------------------------------------------------------------------------ loading
def load(chip, run):
    d = os.path.join(DATA, chip)
    meta = json.load(open(os.path.join(d, 'run%d_meta.json' % run)))
    job = json.load(open(os.path.join(d, 'run%d_job.json' % run)))
    if meta.get('complete') is not True:
        raise SystemExit('run %d on %s not complete' % (run, chip))
    z = np.load(os.path.join(d, 'run%d_results.npz' % run))
    A = {k.split('__')[0]: np.asarray(z[k]).astype(complex) for k in z.files}
    return meta, job, A


def frame(A, off='OFF|p0', on='ON|p0'):
    c0 = A[off].mean(0); d = A[on].mean(0) - c0
    sep = np.abs(d) / np.sqrt((np.abs(A[off] - c0) ** 2).mean(0) / 2)
    return (lambda x: np.real((x - c0) * np.conj(d)) / np.abs(d) ** 2), sep


def resample(A):
    return {k: v[RNG.integers(0, v.shape[0], v.shape[0])] for k, v in A.items()}


def book(m0, q, nseg):
    on = np.asarray(m0, float) * 1.0; off = 1 - on
    for _ in range(nseg):
        on, off = off + on * (1 - q), on * q
    return on


RW = [30, 100, 300, 600]
SPACINGS = [5, 10, 20, 30, 45, 70, 100, 150]


def surv_curve(lv, copy, t):
    ts = np.array([0.0] + RW)
    ys = [np.ones_like(lv['ON|p0'])] + [np.clip(lv['RULER|ON|n%d|%s' % (w, copy)], 0.02, 1.0) for w in RW]
    L = np.log(np.vstack(ys)); i = int(np.searchsorted(ts, t) - 1); i = int(min(max(i, 0), len(ts) - 2))
    f = (t - ts[i]) / (ts[i + 1] - ts[i])
    return np.exp(L[i] + f * (L[i + 1] - L[i]))


# ------------------------------------------------------------------------------------------------------------------ 4 October sweep (runs 301-303; 361)
def sweep_stats(A):
    """per-read numbers of the sweep: for start OFF/ON, spacing, bookkeeping ONE/CURVE: mean over ruler copies of the mean over readable qubits."""
    al, sep = frame(A); g = sep > 4
    keys = ['OFF|p0', 'ON|p0', 'OFF|p10', 'ON|p10'] + ['%s|s10@%d' % (s, t) for s in ('OFF', 'ON') for t in SPACINGS] \
        + ['RULER|ON|n%d|%s' % (w, c) for w in RW for c in 'abc']
    lv = {k: al(A[k]).mean(0) for k in keys}
    M = lambda x: float(np.nanmean(x[g]))
    out = {'_bits': int(g.sum())}
    for s in ('OFF', 'ON'):
        out[(s, 'p10')] = M(lv['%s|p10' % s] - lv['%s|p0' % s])
        start = np.zeros_like(lv['OFF|p0']) if s == 'OFF' else np.ones_like(lv['OFF|p0'])
        for t in SPACINGS:
            v = lv['%s|s10@%d' % (s, t)]
            for c in 'abc':
                q1 = np.clip(lv['RULER|ON|n300|%s' % c], 0.02, 0.98) ** (t / 300.0)
                out[(s, t, 'ONE', c)] = M(v - book(start, q1, 10))
                out[(s, t, 'CURVE', c)] = M(v - book(start, surv_curve(lv, c, t), 10))
    return out


def sweep_chip(chip, runs):
    reads = []
    for r in runs:
        meta, job, A = load(chip, r)
        s = sweep_stats(A)
        B = [sweep_stats(resample(A)) for _ in range(NB)]
        w = {k: float(np.nanstd([b[k] for b in B])) for k in s if k != '_bits'}
        reads.append((r, job, s, w))
    return reads


def comb(reads, keys):
    """sumtime.py rule: value = mean over reads of (mean over keys); width = sqrt(sum of reads' key-averaged widths^2)/n."""
    v = [np.mean([s[k] for k in keys]) for _, _, s, _ in reads]
    wd = [np.mean([w[k] for k in keys]) for _, _, s, w in reads]
    return float(np.mean(v)), float(math.sqrt(np.sum(np.square(wd))) / len(reads))


def sw(reads, s, t, book_):
    return comb(reads, [(s, t, book_, c) for c in 'abc'])


def sw_p10(reads, s):
    return comb(reads, [(s, 'p10')])


# ------------------------------------------------------------------------------------------------------------------ 2 October push runs (91-93; 141)
def push_stats(A):
    al, sep = frame(A); g = sep > 4
    lv = {k: al(v).mean(0) for k, v in A.items()}
    q = np.clip(lv['ON|n300'], 0.02, 0.98) ** 0.1
    out = {'_bits': int(g.sum())}
    for s in ('OFF', 'ON', 'MID'):
        base = lv['%s|p0' % s]
        start = np.zeros_like(base) if s == 'OFF' else (np.ones_like(base) if s == 'ON' else base)
        for n in ('p2', 'p10'):
            out['%s %s' % (s, n)] = float((lv['%s|%s' % (s, n)] - base)[g].mean())
        out['%s s10' % s] = float((lv['%s|s10' % s] - book(start, q, 10))[g].mean())
    return out


def push_read(chip, run):
    meta, job, A = load(chip, run)
    s = push_stats(A); B = [push_stats(resample(A)) for _ in range(NB)]
    w = {k: float(np.std([b[k] for b in B])) for k in s if k != '_bits'}
    return job, s, w


def three(vals, widths):
    v = np.array(vals); w = np.array(widths)
    return float(v.mean()), float(math.sqrt((w ** 2).sum()) / len(v))


# ------------------------------------------------------------------------------------------------------------------ run 121 (2 Oct, no-name shot)
def noname_stats(A):
    al, sep = frame(A, 'F0', 'F1'); g = sep > 4
    lv = {k: al(v).mean(0) for k, v in A.items()}
    R = {L: lv['R%d' % L] for L in (10, 30, 90, 300)}
    one = lambda L: np.clip(lv['R300'], 0.02, 0.98) ** (L / 300.0)
    out = {'_bits': int(g.sum())}
    for L in (10, 30, 90, 300):
        out['ON %d' % L] = float(R[L][g].mean())
    for L in (10, 30, 90):
        out['dev %d' % L] = float((R[L] - one(L))[g].mean())
    # A30_30: start OFF, ten intervals of 30 us, an X after each (stretches [30,30]*5 in noname.specs: OFF-stretch 30, ON-stretch 30)
    out['A30_30 level'] = float(lv['A30_30'][g].mean())
    out['A30_30 minus one-clock'] = float((lv['A30_30'] - book(np.zeros_like(lv['F0']), one(30), 10))[g].mean())
    return out


# ------------------------------------------------------------------------------------------------------------------ gate fit (runs 91-93)
def gate_fit(runs):
    acc = {}
    for r in runs:
        meta, job, A = load('ibm_pittsburgh', r)
        al, sep = frame(A); g = np.where(sep > 4)[0]
        lv = {k: al(v) for k, v in A.items()}
        for arm, n in (('OFF', 2), ('OFF', 10), ('ON', 2), ('ON', 10)):
            a = lv['%s|p%d' % (arm, n)]; b = lv['%s|p0' % arm]
            d = (a.mean(0) - b.mean(0))[g]
            se2 = (a.std(0, ddof=1) ** 2 / a.shape[0] + b.std(0, ddof=1) ** 2 / b.shape[0])[g]
            acc.setdefault((arm, n), []).append((float(d.mean()), float(math.sqrt(se2.sum()) / len(g))))
    val = {k: (float(np.mean([x for x, _ in l])), float(math.sqrt(sum(w * w for _, w in l)) / len(l))) for k, l in acc.items()}
    pts = [(2, val[('OFF', 2)]), (2, (-val[('ON', 2)][0], val[('ON', 2)][1])), (10, val[('OFF', 10)]), (10, (-val[('ON', 10)][0], val[('ON', 10)][1]))]
    A11 = A12 = A22 = B1 = B2 = 0.0
    for N, (y, s) in pts:
        w = 1 / s ** 2; f1 = N * N / 4.0; f2 = float(N)
        A11 += w * f1 * f1; A12 += w * f1 * f2; A22 += w * f2 * f2; B1 += w * f1 * y; B2 += w * f2 * y
    det = A11 * A22 - A12 * A12
    e2 = (B1 * A22 - B2 * A12) / det; b = (A11 * B2 - A12 * B1) / det
    se_e2 = math.sqrt(A22 / det); se_b = math.sqrt(A11 / det)
    eps = math.sqrt(e2)
    return val, eps, se_e2 / (2 * eps), b, se_b


def equator_loss(N, eps):
    """ideal two-level model: perfect sqrt(X) to the equator, then N rotations about x by (pi - eps), no decay between. Returns P1(after) - P1(before)."""
    def rx(th):
        c, s = math.cos(th / 2), math.sin(th / 2); return np.array([[c, -1j * s], [-1j * s, c]])
    psi = rx(math.pi / 2) @ np.array([1, 0], complex)
    p_before = abs(psi[1]) ** 2
    for _ in range(N):
        psi = rx(math.pi - eps) @ psi
    return abs(psi[1]) ** 2 - p_before


# ------------------------------------------------------------------------------------------------------------------ report plumbing
ROWS = []


def row(where, what, paper, got, dec, kind='value', note=''):
    """paper/got in the same units; dec = decimals printed in the paper. status: ok if round(got, dec) == paper (value) or within 15% (width)."""
    if kind == 'width':
        st = 'ok(boot)' if abs(got - paper) <= 0.15 * paper + 0.5 * 10 ** -dec * 0.0 else 'DIFFERS'
        if st == 'DIFFERS' and abs(got - paper) <= 0.5 * 10 ** -dec:
            st = 'ok(boot)'
    else:
        st = 'ok' if abs(round(got + 0.0, dec) - paper) < 0.5 * 10 ** -(dec + 3) else 'DIFFERS'
    ROWS.append((where, what, paper, got, dec, st, note))


def show_rows():
    P('%-9s %-58s %12s %12s  %-9s %s' % ('where', 'quantity', 'paper', 'recomputed', 'status', 'note'))
    for where, what, paper, got, dec, st, note in ROWS:
        P('%-9s %-58s %12s %12s  %-9s %s' % (where, what, ('%+.*f' % (dec, paper)) if paper is not None else '-', ('%+.*f' % (dec + 2, got)) if got is not None else '-', st, note))


def main():
    P('verify_xtrain.py -- Sec. sec:xtrain numbers recomputed from raw IBM job results (data/). Bootstrap NB=%d, seed 20261022.' % NB)
    P('')
    P('JOBS USED (chip run job_id submitted):')
    jobs = {}
    for chip, runs in (('ibm_pittsburgh', (91, 92, 93, 121, 361)), ('ibm_marrakesh', (141, 301, 302, 303)), ('ibm_kingston', (141, 301, 302, 303))):
        for r in runs:
            j = json.load(open(os.path.join(DATA, chip, 'run%d_job.json' % r)))
            jobs[(chip, r)] = j
            P('  %-15s run %-4d %s  %s' % (chip, r, j['job_id'], j['submitted_at']))
    P('')

    # ---- 2 October push runs
    pit = [push_read('ibm_pittsburgh', r) for r in (91, 92, 93)]
    marr141 = push_read('ibm_marrakesh', 141)
    king141 = push_read('ibm_kingston', 141)
    pm = lambda key: three([p[1][key] for p in pit], [p[2][key] for p in pit])
    # ---- run 121
    meta121, job121, A121 = load('ibm_pittsburgh', 121)
    s121 = noname_stats(A121); B121 = [noname_stats(resample(A121)) for _ in range(NB)]
    w121 = {k: float(np.std([b[k] for b in B121])) for k in s121 if k != '_bits'}

    # ---- line 229 (Shape paragraph)
    for L, pv in ((10, 0.9463), (30, 0.8609), (90, 0.6846), (300, 0.3286)):
        row('l.229', 'ON left alone %d us level (run 121)' % L, pv, s121['ON %d' % L], 4)
    for L, pv in ((10, -0.0141), (30, -0.0259), (90, -0.0185)):
        row('l.229', 'ON %d us minus exponential through 300 us' % L, pv, s121['dev %d' % L], 4)
    ROWS.append(('l.229', 'readable qubits in run 121 ("150 qubits")', 150, float(s121['_bits']), 0, 'ok' if s121['_bits'] == 150 else 'DIFFERS', 'of %d' % A121['F0'].shape[1]))

    # ---- sec:xtrain para 1, 2 October spaced trains
    m, w = pm('OFF s10'); row('xtrain', 'pittsburgh 3 jobs, start 0, 10 gates/300 us, vs ONE', 0.0200, m, 4); row('xtrain', '   width', 0.0016, w, 4, 'width')
    pitt_mean_off, pitt_w_off = m, w
    row('xtrain', 'pittsburgh 4th job (run 121), start 0, vs ONE', 0.0126, s121['A30_30 minus one-clock'], 4)
    row('xtrain', '   width', 0.0032, w121['A30_30 minus one-clock'], 4, 'width')
    nsig = abs(pitt_mean_off - s121['A30_30 minus one-clock']) / math.hypot(pitt_w_off, w121['A30_30 minus one-clock'])
    ROWS.append(('xtrain', '"about two widths apart" (diff / combined width)', 2.0, nsig, 0, 'ok' if round(nsig) == 2 else 'DIFFERS', 'ratio, not a level'))
    j, s, w = marr141; row('xtrain', 'marrakesh (1 job), start 0, vs ONE', 0.0274, s['OFF s10'], 4); row('xtrain', '   width', 0.0028, w['OFF s10'], 4, 'width')
    marr_off, marr_w = s['OFF s10'], w['OFF s10']
    j, s, w = king141; row('xtrain', 'kingston (1 job), start 0, vs ONE', 0.0287, s['OFF s10'], 4); row('xtrain', '   width', 0.0028, w['OFF s10'], 4, 'width')
    m, w = pm('ON s10'); row('xtrain', 'pittsburgh start 1, vs ONE', 0.0005, m, 4); row('xtrain', '   width', 0.0013, w, 4, 'width')
    j, s, w = marr141; row('xtrain', 'marrakesh start 1', -0.0016, s['ON s10'], 4); row('xtrain', '   width', 0.0021, w['ON s10'], 4, 'width')
    j, s, w = king141; row('xtrain', 'kingston start 1', -0.0124, s['ON s10'], 4); row('xtrain', '   width', 0.0024, w['ON s10'], 4, 'width')
    d = marr_off - pitt_mean_off; dw = math.hypot(marr_w, pitt_w_off)
    row('xtrain', 'marrakesh minus pittsburgh excess', 0.0074, d, 4); row('xtrain', '   width', 0.0032, dw, 4, 'width')

    # ---- back-to-back, 2 October
    m, w = pm('OFF p10'); row('xtrain', 'pittsburgh 2 Oct back-to-back ten gates, start 0', 0.023, m, 3, note='paper: "+-0.023"')
    m, w = pm('ON p10'); row('xtrain', 'pittsburgh 2 Oct back-to-back ten gates, start 1', -0.023, m, 3, note='paper: "+-0.023"')
    j, s, w = marr141
    row('xtrain', 'marrakesh 2 Oct back-to-back ten gates, start 1', -0.0148, s['ON p10'], 4, note='paper "+-0.0148"; start 0 reads %+.4f' % s['OFF p10'])

    # ---- 4 October sweep
    marr = sweep_chip('ibm_marrakesh', (301, 302, 303)); king = sweep_chip('ibm_kingston', (301, 302, 303)); pitsw = sweep_chip('ibm_pittsburgh', (361,))
    CH = (('marrakesh', marr), ('kingston', king), ('pittsburgh', pitsw))
    P('\n4 October readable qubits per job: ' + '; '.join('%s %s' % (n, [r[2]['_bits'] for r in rs]) for n, rs in CH))
    for (n, rs), pv, pw in zip(CH, (0.0302, 0.0325, 0.0171), (0.0020, 0.0019, 0.0028)):
        m, w = sw(rs, 'OFF', 30, 'ONE'); row('xtrain', '%s start 0 excess vs ONE at 30 us' % n, pv, m, 4); row('xtrain', '   width', pw, w, 4, 'width')
    for (n, rs), pv in zip(CH, (0.0575, 0.0576, 0.0267)):
        m, w = sw(rs, 'OFF', 5, 'ONE'); row('xtrain', '%s start 0 excess vs ONE at 5 us' % n, pv, m, 4)
    for (n, rs), pv in zip(CH, (-0.0003, 0.0058, 0.0050)):
        m, w = sw(rs, 'OFF', 150, 'ONE'); row('xtrain', '%s start 0 excess vs ONE at 150 us' % n, pv, m, 4)
    # start 1 ratio
    ratios = {}
    for n, rs in CH:
        for t in (5, 10, 20):
            ratios[(n, t, 'ONE')] = -sw(rs, 'ON', t, 'ONE')[0] / sw(rs, 'OFF', t, 'ONE')[0]
            ratios[(n, t, 'CURVE')] = -sw(rs, 'ON', t, 'CURVE')[0] / sw(rs, 'OFF', t, 'CURVE')[0]
    vals = [ratios[(n, t, 'ONE')] for n, _ in CH for t in (5, 10, 20)]
    ROWS.append(('xtrain', 'start-1 / start-0 excess at 5,10,20 us vs ONE: min', 0.35, min(vals), 2, 'ok' if round(min(vals), 2) == 0.35 else 'DIFFERS', 'paper: "35 to 100%"'))
    ROWS.append(('xtrain', 'start-1 / start-0 excess at 5,10,20 us vs ONE: max', 1.00, max(vals), 2, 'ok' if round(max(vals), 2) == 1.00 else 'DIFFERS', 'paper: "35 to 100%"'))
    P('\nstart-1 over start-0 excess (minus sign dropped), against ONE (the only bookkeeping for which the ratio is defined; against CURVE the start-0 excess crosses zero):')
    for n, _ in CH:
        P('  %-11s 5 us %.3f   10 us %.3f   20 us %.3f' % ((n,) + tuple(ratios[(n, t, 'ONE')] for t in (5, 10, 20))))
    # back-to-back 4 Oct
    offs = [sw_p10(rs, 'OFF')[0] for _, rs in CH]; ons = [sw_p10(rs, 'ON')[0] for _, rs in CH]
    lo, hi = min(offs) * 100, max(offs) * 100
    ROWS.append(('xtrain', '4 Oct ten back-to-back, start 0: min %', 1.0, lo, 1, 'ok' if round(lo + 1e-9, 1) == 1.0 else 'DIFFERS', 'paper: "+1.0 to +1.2%"'))
    ROWS.append(('xtrain', '4 Oct ten back-to-back, start 0: max %', 1.2, hi, 1, 'ok' if round(hi + 1e-9, 1) == 1.2 else 'DIFFERS', 'paper: "+1.0 to +1.2%"'))
    lo, hi = min(ons) * 100, max(ons) * 100
    ROWS.append(('xtrain', '4 Oct ten back-to-back, start 1: |max| %', 1.1, -lo, 1, 'ok' if round(-lo + 1e-9, 1) == 1.1 else 'DIFFERS', 'paper: "-1.0 to -1.1%"'))
    ROWS.append(('xtrain', '4 Oct ten back-to-back, start 1: |min| %', 1.0, -hi, 1, 'ok' if round(-hi + 1e-9, 1) == 1.0 else 'DIFFERS', 'paper: "-1.0 to -1.1%"'))
    row('xtrain', 'pittsburgh 4 Oct back-to-back ten gates, start 0', 0.010, sw_p10(pitsw, 'OFF')[0], 3, note='paper: "+-0.010"')
    row('xtrain', 'pittsburgh 4 Oct back-to-back ten gates, start 1', -0.010, sw_p10(pitsw, 'ON')[0], 3, note='paper: "+-0.010"')

    # ---- against the job's own decay curve
    def span(rs, s, ts):
        v = [sw(rs, s, t, 'CURVE')[0] for t in ts]; return min(v), max(v), v
    T20 = (20, 30, 45, 70, 100, 150); T30 = (30, 45, 70, 100, 150)
    for n, rs in CH[1:]:
        for s in ('OFF', 'ON'):
            lo, hi, v = span(rs, s, T20)
            ROWS.append(('xtrain', '%s start %s vs CURVE, 20-150 us: max |.|' % (n, '0' if s == 'OFF' else '1'), 0.007, max(abs(lo), abs(hi)), 3, 'ok' if max(abs(lo), abs(hi)) <= 0.0070 + 1e-9 else 'DIFFERS', 'paper: "within +-0.007"; range %+.4f..%+.4f' % (lo, hi)))
            lo, hi, v = span(rs, s, T30)
            ROWS.append(('xtrain', '%s start %s vs CURVE, 30-150 us: max |.|' % (n, '0' if s == 'OFF' else '1'), 0.007, max(abs(lo), abs(hi)), 3, 'ok' if max(abs(lo), abs(hi)) <= 0.0070 + 1e-9 else 'DIFFERS', 'paper: "within +-0.007" at 30 us and longer; range %+.4f..%+.4f' % (lo, hi)))
    lo, hi, v = span(marr, 'OFF', T20)
    ROWS.append(('xtrain', 'marrakesh start 0 vs CURVE, 20-150 us: least negative', -0.012, hi, 3, 'ok' if round(hi, 3) == -0.012 else 'DIFFERS', 'paper: "-0.012 to -0.018"'))
    ROWS.append(('xtrain', 'marrakesh start 0 vs CURVE, 20-150 us: most negative', -0.018, lo, 3, 'ok' if round(lo, 3) == -0.018 else 'DIFFERS', 'paper: "-0.012 to -0.018"'))
    lo30, hi30, _ = span(marr, 'OFF', T30)
    ROWS.append(('xtrain', 'marrakesh start 0 vs CURVE, 30-150 only: least negative', -0.012, hi30, 3, 'info', 'info: the -0.012 end needs the 20 us (interpolated) point'))
    loB, hiB, vB = span(marr, 'ON', T20)
    ROWS.append(('xtrain', 'marrakesh start 1 vs CURVE, 20-150 us: most negative', -0.016, loB, 3, 'ok' if round(loB, 3) == -0.016 else 'DIFFERS', 'paper: "by up to -0.016"; max positive %+.4f' % hiB))
    ROWS.append(('abstract', 'marrakesh over-correction 20-150 us, start 0 least/most (%)', None, None, 0, 'info', '%.2f to %.2f %%; start 1 up to %.2f %% (paper: 1 to 2%%)' % (-hi * 100, -lo * 100, -loB * 100)))
    # residuals at 5 and 10 us
    for n, rs, p5, p10_ in (('marrakesh', marr, 0.0191, 0.0040), ('kingston', king, 0.0337, 0.0195), ('pittsburgh', pitsw, 0.0108, 0.0073)):
        row('xtrain', '%s start 0 vs CURVE at 5 us' % n, p5, sw(rs, 'OFF', 5, 'CURVE')[0], 4)
        row('xtrain', '%s start 0 vs CURVE at 10 us' % n, p10_, sw(rs, 'OFF', 10, 'CURVE')[0], 4)
        sg = [sw(rs, 'ON', t, 'CURVE')[0] for t in (5, 10)]
        ROWS.append(('xtrain', '%s start 1 vs CURVE at 5, 10 us (paper: opposite sign)' % n, None, sg[0], 4, 'ok' if all(x < 0 for x in sg) else 'DIFFERS', '5 us %+.4f, 10 us %+.4f' % (sg[0], sg[1])))

    # ---- gate fit and equator prediction
    val, eps, seps, b, seb = gate_fit((91, 92, 93))
    row('gate', 'ten back-to-back gates, start 0 (%)', 2.29, 100 * val[('OFF', 10)][0], 2)
    row('gate', 'ten back-to-back gates, start 1 (%)', -2.32, 100 * val[('ON', 10)][0], 2)
    row('gate', 'X shortfall (deg)', 1.47, math.degrees(eps), 2); row('gate', '   width (deg)', 0.12, math.degrees(seps), 2, 'width')
    row('gate', '4b (% per gate)', 0.261, 400 * b, 3); row('gate', '   width of 4b (%)', 0.104, 400 * seb, 3, 'width')
    row('gate', 'upper edge 4(b + 2 se) (%)', 0.47, 400 * (b + 2 * seb), 2)
    ph = pm('MID p2'); pt = pm('MID p10')
    row('gate', 'equator start, 2 gates, measured loss', 0.0262, -ph[0], 4)
    row('gate', 'equator start, 10 gates, measured loss', 0.114, -pt[0], 3)
    for N, pv, pw in ((2, 0.026, 0.002), (10, 0.127, 0.012)):
        pred = -equator_loss(N, eps)
        dpred = abs(equator_loss(N, eps + 1e-6) - equator_loss(N, eps - 1e-6)) / 2e-6 * seps
        row('gate', 'shortfall-alone prediction, %d gates' % N, pv, pred, 3)
        row('gate', '   width from the shortfall width', pw, dpred, 3, 'width')
        if N == 10:
            need = pw / (abs(equator_loss(N, eps + 1e-6) - equator_loss(N, eps - 1e-6)) / 2e-6)
            P('    info: the paper width %.3f on the 10-gate prediction needs a shortfall width of %.3f deg (fit gives %.3f deg)' % (pw, math.degrees(need), math.degrees(seps)))
    P('')
    show_rows()
    n_bad = sum(1 for r in ROWS if r[5] == 'DIFFERS')
    P('\n%d rows, %d marked DIFFERS.' % (len(ROWS), n_bad))
    P('\nNOTES')
    P(' 1. Rows marked DIFFERS are listed first-hand above; none was adjusted. The bookkeeping conventions are exactly those of readclose.py / readtime.py / readnoname.py.')
    P(' 2. "start-1 / start-0 ... 35 to 100%": recomputed range is 37 to 97% (pittsburgh 20 us and 5 us are the ends); the paper rounds outward.')
    P(' 3. Ten-gate prediction width: 0.012 in the paper; the fitted shortfall width (%.3f deg) gives %.4f. 0.012 needs a shortfall width of about 0.142 deg;' % (math.degrees(seps), dpred))
    P('    the notebook (section 149, 1.48 deg) printed 0.128 +- 0.012 and 0.0258 +- 0.0024, i.e. a width of about 0.14 deg, not the 0.12 deg now in the paper.')
    P(' 4. Pittsburgh 4 October: runs 362 and 363 (2 more jobs, same circuits, 15:33 and 15:34 EDT) exist in the record; the paper uses run 361 only ("one on pittsburgh"), as here.')
    P(' 5. Not numeric, not recomputed: "the population in the band where a |2> level at mid-line would read is flat across arms" and "the bookkeeping removes most of the excess".')
    P(' 6. Widths are bootstrap values and differ from the reference scripts by their own resampling noise (up to ~12% here); values (means) are not resampled and match exactly.')
    open(os.path.join(HERE, 'results_xtrain.text'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
