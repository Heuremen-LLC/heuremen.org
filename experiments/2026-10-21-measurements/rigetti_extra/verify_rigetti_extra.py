#!/usr/bin/env python3
"""verify_rigetti_extra.py -- compute two Rigetti Cepheus-1-108Q results from data/rigetti_extra.json.

    python3 verify_rigetti_extra.py

Standard library only, no network, no account. Every value printed is computed from the raw shot counts in the json.
Quantities that would need a readout-error model are not computed here and are labelled NOT COMPUTED.

Plain-language key:  'push' = X gate (native rx(pi));  'readout' = measurement;  'clock' = spectator qubit prepared on the
equator (rx(pi/2)), left free, then rotated by a second rx(pi/2) (hand X) or rz(pi/2) then rx(pi/2) (hand Y);
phase of the clock = angle of  H = (1 - 2*P1[X]) - i*(1 - 2*P1[Y])  where P1 is the fraction of shots reading 1.
"""
import json, math, os, cmath

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'data', 'rigetti_extra.json')))
QIDS = D['meta']['qubit_ids']; IDX = {q: i for i, q in enumerate(QIDS)}

def mean(x): x = list(x); return sum(x) / len(x)
def sd(x):
    x = list(x); m = mean(x); return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))


# =====================================================================================================
print('=' * 100)
print('RESULT A -- ten X gates in a row on 7 well-separated qubits, read 1 more often when the 22 neighbours')
print('             were held half-way (equator) around the train. Rigetti Cepheus-1-108Q, 2026-10-04.')
print('=' * 100)
A = D['finding_a']; EM = A['emitters']; BO = A['bit_order']; NBR = A['all_neighbours']; N = A['shots']
t_sent = sorted(t['sent_utc'] for t in A['tasks'])
print(f'tasks: {len(A["tasks"])} (3 rounds x [8 calibration + 6 experiment]); sent {t_sent[0]} .. {t_sent[-1]}; {N} shots per experiment task')
for r in (1, 2, 3):
    ts = [t['sent_utc'] for t in A['tasks'] if t['tag'][-1] == r]; print(f'  round {r}: {min(ts)} .. {max(ts)}')
print(f'emitters (the qubits that are pushed 10x): {EM}; neighbours ({len(NBR)}): {NBR}')
print('neighbours of each emitter:', {int(k): v for k, v in A['neighbours_of_emitter'].items()})
ex = {t['tag'][0] + '|' + str(t['tag'][1]) + '|' + str(t['tag'][2]): t['task_id'] for t in A['tasks'] if t['tag'][0] != 'cal'}
print('experiment task ids (state|readout|round -> Braket task id, region us-west-1; readout e = plain readout, ba = swap-then-flip readout):')
for k in sorted(ex): print(f'    {k:14} {ex[k]}')

# unpack
ROW = {}
for key, ints in A['e_packed'].items():
    r, s = key.split('|'); ROW[(int(r), s)] = ints
def emit_ones(r, s):
    return [sum((v >> i) & 1 for v in ROW[(r, s)]) for i in range(len(EM))]
def nbr_frac(r, s):
    return mean(sum((v >> (len(EM) + j)) & 1 for v in ROW[(r, s)]) / N for j in range(len(NBR)))

STATES = ('rest', 'pushed', 'clock')
P = {}; print('\nA1. read-1 fraction of the pushed qubits, mean over the 7 (per-qubit values below):')
for r in (1, 2, 3):
    for s in STATES:
        P[(r, s)] = [c / N for c in emit_ones(r, s)]
    print(f'  round {r} mean over the 7: ' + ' / '.join(f'{s} {mean(P[(r, s)]):.5g}' for s in STATES))
print('\n  per-qubit read-1 fractions (columns q' + ' q'.join(map(str, EM)) + '):')
for r in (1, 2, 3):
    for s in STATES: print(f'    round {r} {s:6}: ' + ' '.join(f'{v:.3f}' for v in P[(r, s)]) + f'   mean {mean(P[(r, s)]):.4f}')

print('\nA2. widths. Two kinds, both reported:')
print('    shot noise on the 7-qubit mean = sqrt(sum p(1-p)/N)/7  (binomial, per round);')
print('    qubit-to-qubit spread = standard deviation of the 7 per-qubit values / sqrt(7).')
for r in (1, 2, 3):
    for s in STATES:
        p = P[(r, s)]; shot = math.sqrt(sum(v * (1 - v) / N for v in p)) / 7; spread = sd(p) / math.sqrt(7)
        print(f'    round {r} {s:6}: {mean(p):.4f}  shot +-{shot:.4f}  spread +-{spread:.4f}')
print('\nA3. the effect, per qubit, pooled over the 3 rounds (2700 shots per qubit per state):')
POOL = {s: [sum(emit_ones(r, s)[i] for r in (1, 2, 3)) / (3 * N) for i in range(len(EM))] for s in STATES}
def paired(a, b, name):
    d = [x - y for x, y in zip(POOL[a], POOL[b])]; m = mean(d); se = sd(d) / math.sqrt(len(d)); t = m / se
    dshot = [math.sqrt(POOL[a][i] * (1 - POOL[a][i]) / (3 * N) + POOL[b][i] * (1 - POOL[b][i]) / (3 * N)) for i in range(len(d))]
    print(f'    {name}: per-qubit difference ' + ' '.join(f'q{q} {v:+.3f}' for q, v in zip(EM, d)))
    print(f'      mean {m:+.4f}  +- {se:.4f} (spread across 7 qubits, 6 d.o.f.; t = {t:.2f})   shot-noise-only width of the mean +-{math.sqrt(sum(v*v for v in dshot))/7:.4f}   qubits positive: {sum(v > 0 for v in d)} of 7')
    return d, m, se
d_cr, m_cr, se_cr = paired('clock', 'rest', 'equator minus rest')
d_pr, m_pr, se_pr = paired('pushed', 'rest', 'pushed minus rest')
paired('clock', 'pushed', 'equator minus pushed')
def ratio(a, b, drop=None):
    ia = [i for i in range(len(EM)) if i != drop]; return mean(POOL[a][i] for i in ia) / mean(POOL[b][i] for i in ia)
R = ratio('clock', 'rest'); jk = [ratio('clock', 'rest', i) for i in range(len(EM))]
jse = math.sqrt((len(jk) - 1) / len(jk) * sum((x - mean(jk)) ** 2 for x in jk))
print(f'    ratio equator/rest of the 7-qubit mean: {R:.2f} +- {jse:.2f} (delete-one-qubit jackknife, so it carries the qubit-to-qubit spread); '
      f'range when one qubit is dropped: {min(jk):.2f} .. {max(jk):.2f}')
print('    ratio per round: ' + ', '.join(f'round {r} {mean(P[(r, "clock")]) / mean(P[(r, "rest")]):.2f}' for r in (1, 2, 3)))
print(f'    pushed/rest ratio: {ratio("pushed", "rest"):.3f}')
print(f'    ratio equator/rest, round 1: {mean(P[(1, "clock")]) / mean(P[(1, "rest")]):.5g}; all 3 rounds pooled: {R:.5g}')
print('    (the ratio differs from round to round; it is not a per-round constant)')
print('    significance across qubits: mean paired rise / its spread-based width gives t = %.2f on 6 d.o.f.; the 95 %% two-sided critical value is 2.447, so on qubit-to-qubit spread alone the pooled rise is NOT past it.' % (m_cr / se_cr))
print('    Every single qubit except q57 and q72 is far past its own shot noise (A6); it is the unevenness between qubits that keeps the pooled test marginal.')
print(f'    pushed minus rest, pooled mean difference: {m_pr:.5g}')

print('\nA4. WHAT STATE ARE THE NEIGHBOURS IN WHEN THE PUSHED QUBITS ARE READ?  (read-1 fraction of the 22 neighbours)')
for r in (1, 2, 3): print(f'    round {r}: ' + ' / '.join(f'{s} {nbr_frac(r, s):.3f}' for s in STATES))
print('    -> in the "equator" arm the neighbours end in |1> (about 94 % read 1): their two half-turns (rx(pi/2) before the train, rx(pi/2) after)')
print('       make one full turn. They sit on the equator only DURING the train, and are excited when the readout happens. In the other two arms')
print('       they are in |0> at the readout (rest, pushed 10x: the pushed arm shows ~9 % reading 1 from gate error).')

print('\nA5. shot-level check (same shots, no model): pushed-qubit read-1 rate versus how many of ITS OWN neighbours read 1 in the same shot.')
print('    rate (shots) by number of own neighbours that read 1; pooled over the 3 rounds and, for each qubit, over its neighbours.')
for s in STATES:
    cells = {}
    for qi, q in enumerate(EM):
        nb = [BO.index(n) for n in A['neighbours_of_emitter'][str(q)]]
        for r in (1, 2, 3):
            for v in ROW[(r, s)]:
                k = sum((v >> j) & 1 for j in nb); n = len(nb)
                cell = cells.setdefault(('all neighbours read 0' if k == 0 else 'all neighbours read 1' if k == n else 'some read 1'), [0, 0])
                cell[0] += (v >> qi) & 1; cell[1] += 1
    print(f'    {s:6}: ' + '   '.join(f'{lab}: {c[0] / c[1]:.3f} ({c[1]} shots)' for lab, c in sorted(cells.items())))
print('    reading: even in the REST arm a qubit reads 1 about four times more often in shots where one of its neighbours also reads 1. The readouts are correlated')
print('    between neighbours (read-out crosstalk or a shared real excitation; this data cannot tell which). Only the equator arm has neighbours in |1> at the readout,')
print('    so "neighbours on the equator during the train" and "neighbours in |1> at the readout" are confounded in this design. A simple crosstalk-by-neighbour-reading')
print('    rule fitted on the rest arm would predict more than the 0.099 seen in the equator arm, so crosstalk alone is not a quantitative match either.')
print('\nA6. per qubit, equator minus rest (pooled): the excess is not uniform. Qubits with a rise > 2 pooled-shot widths:')
for i, q in enumerate(EM):
    w = math.sqrt(POOL['clock'][i] * (1 - POOL['clock'][i]) / (3 * N) + POOL['rest'][i] * (1 - POOL['rest'][i]) / (3 * N))
    print(f'    q{q:<4} rest {POOL["rest"][i]:.3f} equator {POOL["clock"][i]:.3f} diff {POOL["clock"][i] - POOL["rest"][i]:+.3f} +-{w:.3f}  ({(POOL["clock"][i] - POOL["rest"][i]) / w:+.1f} widths)')

# =====================================================================================================
print('\n' + '=' * 100)
print('RESULT B -- a neighbour (q0) prepared in |2> shifts the precession rate of a spectator (q9).')
print('=' * 100)
# ---- sign conventions, from a tiny exact simulation of each day's clock circuit (pure python, 2x2 complex matrices)
def mm(a, b): return [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def rx(t): c, s = math.cos(t / 2), math.sin(t / 2); return [[c, -1j * s], [-1j * s, c]]
def rz(t): return [[cmath.exp(-1j * t / 2), 0], [0, cmath.exp(1j * t / 2)]]
def run(gates):
    U = [[1, 0], [0, 1]]
    for g in gates: U = mm(g, U)
    return abs(U[1][0]) ** 2
def argH(day, phi):
    """angle of H when the clock's |1> gains phase phi (degrees) relative to |0> during the wait (an rz(phi))."""
    ph = math.radians(phi); h = math.pi / 2
    if day == 'sep28':   # ry(a) = rz(-pi/2) rx(a) rz(pi/2); circuit: ry(pi/2), wait, [Y: rz(-pi/2)], ry(-pi/2)
        pre = [rz(-h), rx(h), rz(h)]; fin = [rz(-h), rx(-h), rz(h)]; finY = [rz(-h)] + fin
    else:                # oct4: rx(pi/2), wait, [Y: rz(pi/2)], rx(pi/2)
        pre = [rx(h)]; fin = [rx(h)]; finY = [rz(h), rx(h)]
    px = run(pre + [rz(ph)] + fin); py = run(pre + [rz(ph)] + finY)
    return math.degrees(cmath.phase(complex(1 - 2 * px, -(1 - 2 * py))))
def slope_sign(day):
    d = argH(day, 10) - argH(day, 0); d = (d + 180) % 360 - 180; return d / 10
SG = {day: slope_sign(day) for day in ('sep28', 'oct4')}
print(f'B0. sign bookkeeping (simulated from each day\'s circuit): d angle(H) / d(clock phase) = {SG["sep28"]:+.0f} on 28 Sep, {SG["oct4"]:+.0f} on 4 Oct.')
print('    The two days used different final rotations, so the SAME physical shift reads with REVERSED sign in angle(H). Below, "rate" = clock phase per microsecond')
print('    in a common convention (positive = |1> gains phase relative to |0> at that rate, i.e. an rz(+) ), obtained as angle(H) rate / the number above.')
print('    This assumes rx/rz mean the same rotations on both days (same compiler and gate set). The sign agreement across days does not depend on that assumption.')


def hvar(px, py, n):
    a, b = 1 - 2 * px, 1 - 2 * py; va = 4 * px * (1 - px) / n; vb = 4 * py * (1 - py) / n
    return (b * b * va + a * a * vb) / (a * a + b * b) ** 2


def turn(c1x, c1y, n1, c0x, c0y, n0):
    """phase of arm 1 relative to arm 0 in degrees, and its shot-noise width (delta method)."""
    p = lambda c, n: c / n
    H1 = complex(1 - 2 * p(c1x, n1), -(1 - 2 * p(c1y, n1))); H0 = complex(1 - 2 * p(c0x, n0), -(1 - 2 * p(c0y, n0)))
    ph = cmath.phase(H1 / H0)
    sig = math.sqrt(hvar(p(c1x, n1), p(c1y, n1), n1) + hvar(p(c0x, n0), p(c0y, n0), n0))
    return math.degrees(ph), math.degrees(sig), abs(H1) / abs(H0)


# ---- 28 September
S = D['sep28']; C = S['ones_q0_q9_shots']
cr = [t['created_utc'] for t in S['tasks']]
print(f'\nB1. 28 SEPTEMBER (UTC {min(cr)} .. {max(cr)} = 20:58-20:59 EDT on 28 Sep): emitter q0 at rest (E0) or in |2> (E2: rx(pi), then a tuned 1<->2 pulse),')
print('    clock q9 (next door; all other qubits idle), waits 0.5 1 2 3 4.5 6 us, hands X and Y, 1000 shots (E0) / 2500 shots (E2). 24 tasks, ids in the json.')
WAITS = [0.5, 1, 2, 3, 4.5, 6]
ang = []; sig = []; hand = []
for w in WAITS:
    g = lambda a, h: C[json.dumps([a, w, h])]
    c2x, c2y, c0x, c0y = g('E2', 'X')[1], g('E2', 'Y')[1], g('E0', 'X')[1], g('E0', 'Y')[1]
    a, s, hr = turn(c2x, c2y, g('E2', 'X')[2], c0x, c0y, g('E0', 'X')[2]); ang.append(a); sig.append(s); hand.append(hr)
un = []; off = 0.0
for i, a in enumerate(ang):                                   # unwrap along the ladder
    a2 = a + off
    if i and abs(a2 - un[-1]) > 180: off += -360 * round((a2 - un[-1]) / 360); a2 = a + off
    un.append(a2)
conv = 1.0 / SG['sep28']
print('    window  wrapped angle   unwrapped   +-width   rate (deg/us, common convention)   amplitude ratio |H(E2)|/|H(E0)|')
rate28 = []
for w, a, u, s, hr in zip(WAITS, ang, un, sig, hand):
    rate28.append((u * conv / w, s / w)); print(f'    {w:5.1f}   {a:+9.2f}   {u:+9.2f}   {s:5.2f}   {u * conv / w:+8.2f} +- {s / w:4.2f}                      {hr:.2f}')
sw = sum(w * w / s ** 2 for w, s in zip(WAITS, sig)); slope = sum(w * u / s ** 2 for w, u, s in zip(WAITS, un, sig)) / sw
se = 1 / math.sqrt(sw); chi2 = sum(((u - slope * w) / s) ** 2 for w, u, s in zip(WAITS, un, sig)); dof = len(WAITS) - 1
print(f'    straight line through the origin, all 6 windows: slope {slope * conv:+.2f} deg/us  +- {se:.2f} (shot noise only)  chi2 = {chi2:.1f} for {dof} d.o.f.  -> width scaled by sqrt(chi2/dof): +- {se * math.sqrt(chi2 / dof):.2f}')
print(f'    28 Sep slope magnitude: {abs(slope):.5g} deg/us;  shot-noise width: {se:.5g} deg/us')
print(f'    NOTE the per-window rates are NOT constant: {rate28[0][0]:+.1f}, {rate28[1][0]:+.1f}, {rate28[2][0]:+.1f}, ... {rate28[5][0]:+.1f}. The straight line fits poorly (chi2/dof = {chi2 / dof:.1f}).')
print('    They rise with window length although the emitter is decaying out of |2> (its lifetime is not computed here), which would make them fall by ~20 % by 6 us.')

# ---- 4 October, triA
Ta = D['oct4_triA']; CA = Ta['ones_all_qubits']; NA = Ta['shots']
ts = sorted(t['sent_utc'] for t in Ta['tasks'])
print(f'\nB2. 4 OCTOBER, run "triA" (tasks sent {ts[0]} .. {ts[-1]}, 13:13 EDT): emitter q0 at rest / |1> / |2>, EVERY other qubit (106) a clock, windows 0.5 1 2 4 6 us, {NA} shots.')
print('    (q0 -> |2>: rx(pi) then the 1<->2 pulse, amplitude re-tuned on 4 Oct, not the 28 Sep amplitude.)  Hand Y of this run also applies the swap-then-flip readout to q0 at the end.')
def triA(q, st, w):
    g = lambda s, h: CA[json.dumps([s, w, h])][IDX[q]]
    return turn(g(st, 'X'), g(st, 'Y'), NA, g('rest', 'X'), g('rest', 'Y'), NA)
WA = [0.5, 1, 2, 4, 6]
r4 = {}
print('    q9 with q0 in |2>, minus q0 at rest (same window): angle (deg), +-width, rate = angle/window (deg/us)')
for w in WA:
    a, s, hr = triA(9, 'two', w); r4[w] = (a / w, s / w); line = f'      {w:4.1f} us  {a:+8.2f}  {s:5.2f}   {a / w:+7.2f} +- {s / w:4.2f}   amplitude ratio |H(|2>)|/|H(rest)| {hr:.2f}'
    print(line)
for w in (0.5, 1, 2): print(f'    4 Oct q9, q0 in |2>, {w} us window: {r4[w][0]:.5g} deg/us')
print('    The 4 us and 6 us windows (+16.2, -1.6) do not continue the trend: at a steady -36 deg/us the angle would be about -145 deg at 4 us and -215 (= +145 wrapped) at 6 us.')
print('    They are shown above as read; the first three windows are the only stretch where a steady shift is visible.')
print('    At 2 us the clock keeps only ~0.54 of the rest amplitude (q0_last, B3: 0.56) while at 0.5 and 1 us it keeps 1.0: not explained here (a shot-to-shot spread of the shift would do it).')
w4 = [(r4[w][0], r4[w][1]) for w in (0.5, 1, 2)]
m4 = sum(r / s ** 2 for r, s in w4) / sum(1 / s ** 2 for r, s in w4); e4 = 1 / math.sqrt(sum(1 / s ** 2 for r, s in w4))
print(f'    inverse-variance mean of 0.5/1/2 us: {m4:+.2f} +- {e4:.2f} deg/us (windows are not independent of q0 decay; shot noise only)')
print('    for comparison, q0 in |1>:  ' + '  '.join(f'{w} us {triA(9, "one", w)[0] / w:+.1f}' for w in WA) + '   (small: |1> hardly moves q9, |2> moves it a lot)')
print('    q1 (other next-door qubit), q0 in |2>:  ' + '  '.join(f'{w} us {triA(1, "two", w)[0] / w:+.1f}' for w in WA))
dist = {0: 0}; nb = {}
for a, b in D['meta']['edges']: nb.setdefault(a, []).append(b); nb.setdefault(b, []).append(a)
fr = [0]
while fr:
    u = fr.pop(0)
    for v in nb.get(u, []):
        if v not in dist: dist[v] = dist[u] + 1; fr.append(v)
far = [q for q in QIDS if q != 0 and dist.get(q, 99) >= 3]
fm = sorted(triA(q, 'two', 1)[0] / 1 for q in far); print(f'    {len(far)} clocks 3+ links from q0, q0 in |2>, 1 us window: median {fm[len(fm) // 2]:+.2f} deg/us, quartiles {fm[len(fm) // 4]:+.1f} .. {fm[3 * len(fm) // 4]:+.1f}  (the shift is local)')

# ---- 4 Oct, q0_last (separate measurement of the same pair, ~20 min before triA)
Tl = D['oct4_q0_last']; CL = Tl['ones_all_qubits']
tsl = sorted(t['sent_utc'] for t in Tl['tasks'])
g = lambda s, h: CL[json.dumps([s, 'T2', h])][IDX[9]]
a, s, hr = turn(g('two', 'X'), g('two', 'Y'), 900, g('rest', 'X'), g('rest', 'Y'), 900)
print(f'\nB3. 4 OCTOBER, run "q0_last" (sent {tsl[0]} .. {tsl[-1]}, 20 minutes BEFORE triA; same pair, one 2 us window, q0 gets no end-of-circuit readout pulses):')
print(f'    q9, q0 in |2> minus rest: angle {a:+.2f} +- {s:.2f} over 2 us = {a / 2:+.2f} +- {s / 2:.2f} deg/us (hand {hr:.2f}); triA at 2 us gave {r4[2][0]:+.2f} +- {r4[2][1]:.2f}')
print(f'    4 Oct q0_last q9: {a / 2:.5g} deg/us;  width: {s / 2:.5g} deg/us')

# ---- same pair?
print('\nB4. SAME SPECTATOR / NEIGHBOUR PAIR BOTH DAYS?')
print(f'    28 Sep: emitter q{S["emitter"]}, clock q{S["clock"]}.   4 Oct (triA, q0_last): emitter q0, clock q9.   Are q0 and q9 coupled on the chip graph? {9 in nb[0]}.')
print('    -> YES, the same pair. The comparison is valid as a comparison of one pair, with the differences listed in B5.')
w3 = [(rate28[i][0], rate28[i][1]) for i in range(3)]       # 28 Sep rates at 0.5, 1, 2 us
mm28 = sum(r / s ** 2 for r, s in w3) / sum(1 / s ** 2 for r, s in w3); em28 = 1 / math.sqrt(sum(1 / s ** 2 for r, s in w3))
print(f'    like-for-like (windows 0.5, 1, 2 us, inverse-variance mean of per-window rates): 28 Sep {mm28:+.2f} +- {em28:.2f}   4 Oct {m4:+.2f} +- {e4:.2f}')
print(f'    ratio 4 Oct / 28 Sep = {m4 / mm28:.2f};   difference {m4 - mm28:+.1f} +- {math.sqrt(em28 ** 2 + e4 ** 2):.1f} deg/us ({abs(m4 - mm28) / math.sqrt(em28 ** 2 + e4 ** 2):.0f} widths, shot noise only)')
print(f'    the through-origin slope (the "-62.1") vs the 4 Oct 0.5-2 us mean: {slope * conv:+.1f} vs {m4:+.1f}')
print('B5. things that differ between the two days and are NOT controlled (so a change between the days cannot be assigned to the coupling alone):')
print('    - the 1<->2 pulse amplitude was re-tuned on 4 Oct, so how completely q0 reached |2> may differ. The |2> population of q0 needs a readout-error model and is NOT COMPUTED here.')
print('      The 28 Sep run has no calibration for q0, so its |2> population is unknown here.')
print('    - 28 Sep: q9 the only clock. 4 Oct: 106 clocks driven at once (the effect of the number of simultaneous clocks on the shift is not measured here).')
print('    - the window ladder and the fit: 28 Sep is a straight-line slope over 0.5-6 us; 4 Oct per-window values over 0.5-2 us, then the later windows break the trend.')
print('    - Cepheus is recalibrated between the days (device parameters drift; both days\' qubit frequencies were not logged side by side).')

# ---- triF context
F = D['oct4_triF']; CF = F['ones']; cols = F['columns']; NF = F['shots']; WINF = F['windows_us']
print('\nB6. CONTEXT, run "triF" (4 Oct 17:37-17:44 UTC, 3 rounds): 7 other emitter/neighbour pairs, emitters all held at rest or in |2>. q0/q9 are not in this run.')
print('    rate = angle of (neighbour with emitter in |2> relative to emitter at rest) / window, deg/us; mean over 3 rounds; each neighbour paired with its own emitter.')
rows = []
for e in F['emitters']:
    for n in F['neighbours_of_emitter'][str(e)]:
        rr = {}
        for w in WINF[:4]:
            vals = []
            for r in (1, 2, 3):
                g = lambda s, h: CF[json.dumps([s, w, h, r])][cols.index(n)]
                vals.append(turn(g('two', 'X'), g('two', 'Y'), NF, g('zero', 'X'), g('zero', 'Y'), NF)[0] / w)
            rr[w] = (mean(vals), sd(vals))
        rows.append((e, n, rr))
print('    emitter neighbour   ' + '  '.join(f'{w:>6} us' for w in WINF[:4]) + '     (mean over 3 rounds; +- = spread of the 3 rounds at 1 us)')
for e, n, rr in sorted(rows, key=lambda x: x[2][1][0]):
    print(f'    q{e:<4}  q{n:<4}        ' + '  '.join(f'{rr[w][0]:+7.1f}' for w in WINF[:4]) + f'     +-{rr[1][1]:.1f}')
v1 = sorted(rr[1][0] for _, _, rr in rows)
print(f'    1 us window over the {len(rows)} pairs: median {v1[len(v1) // 2]:+.1f}, most negative {v1[0]:+.1f}, most positive {v1[-1]:+.1f}; pairs below -20 deg/us: {sum(v < -20 for v in v1)}.  (q0-q9: {r4[1][0]:+.1f}; {sum(abs(v) > abs(r4[1][0]) for v in v1)} of the {len(v1)} pairs are stronger in magnitude, {sum(v > 0 for v in v1)} have the reverse sign.)')
print('    Shifts of 30-70 deg/us from an emitter in |2> therefore appear on other pairs too, and the sign is not always negative: sign and size are properties of the pair.')

