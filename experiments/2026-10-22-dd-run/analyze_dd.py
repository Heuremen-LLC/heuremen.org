"""Analysis of the DD run exactly as PLAN_DD.md. python analyze_dd.py [folder]. Reads job_dd_<chip>.json, raw_dd_<chip>.npz, props_dd_<chip>.json.
The ruler is the deciding runs' (deciding_runs/analyze.py, imported, not copied): every reading a GAP between arms of the same job, per
spectator, median over spectators, 95 % bootstrap (2000 resamples), ABOVE / BELOW / CONSISTENT / CANNOT TELL; copies compared by their
per-spectator difference (excludes 0 -> CANNOT TELL with both copies printed). Signs carried. Chips combined only by inverse variance,
chi2/dof printed. Phase convention of the control run: phase = atan2(<Y>, <X>), gap = wrap(phase_M - phase_D), degrees.
Every reading is printed for all spectators, then WIRED (cz-coupled in props to a measured qubit) and UNWIRED."""
import os, sys, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'deciding_runs'))
import analyze as A
NSTEP = (2, 4, 8, 16)
def dd(J):
    I = J.info; S = np.array(I['S']); E = set(I['E']); n = J.props_n; adj = A.adjacency(J.props, n)
    wired = np.array([len(adj[s] & E) > 0 for s in S]); SUB = (('all', np.ones(len(S), bool)), ('wired', wired), ('unwired', ~wired))
    print(' spectators %d: wired %d, unwired %d | L %s dt, X %s dt, echo %s / %s / %s dt' % (len(S), wired.sum(), (~wired).sum(), I['L_dt'], I['g_dt'], I['d1'], I['d2'], I['d3']))
    def ex(c, d, mode, N, st): return J.mean_pm(c, 'fin', n, dd=d, mode=mode, N=N, setting=st)[S]
    def ph(c, d, mode, N): return np.arctan2(ex(c, d, mode, N, 'YY'), ex(c, d, mode, N, 'XX'))
    def co(c, d, mode, N): return np.hypot(ex(c, d, mode, N, 'YY'), ex(c, d, mode, N, 'XX'))
    # R8 dead-clock guard: coherence (both copies pooled) below 5 x shot noise, shot noise = 1/sqrt(shots per setting, both copies)
    nsh = int(J.rec['shots']) * len(J.idx((0, 1), dd='off', mode='D', N=0, setting='XX'))
    floor = 5 / np.sqrt(nsh); alive = {}
    print(' R8 DEAD-CLOCK GUARD: coherence < 5 x shot noise = %.4f (shots per setting, both copies: %d); dead count of %d (wired / unwired)' % (floor, nsh, len(S)))
    for N in (0,) + NSTEP:
        cells = []
        for d in ('off', 'on'):
            for mode in (('D',) if N == 0 else ('D', 'M')):
                a = co((0, 1), d, mode, N) >= floor; alive[(d, mode, N)] = a
                cells.append('DD %s %s: %d (%d / %d)' % (d, mode, (~a).sum(), (~a & wired).sum(), (~a & ~wired).sum()))
        print('   N=%2d  %s' % (N, ' | '.join(cells)))
    live = {d: np.all([alive[(d, mo, N)] for N in NSTEP for mo in ('D', 'M')], 0) for d in ('off', 'on')}; live['both'] = live['off'] & live['on']
    print('   alive in every arm and N used by the slope: DD off %d, DD on %d, both %d (phase readings use only these)' % (live['off'].sum(), live['on'].sum(), live['both'].sum()))
    def gap(c, d):   # unwrapped along N: phase gap M - D at N = 2, 4, 8, 16, degrees
        g = {N: A.wrap(ph(c, d, 'M', N) - ph(c, d, 'D', N)) for N in NSTEP}; t = {NSTEP[0]: g[NSTEP[0]]}
        for a, b_ in zip(NSTEP[:-1], NSTEP[1:]): t[b_] = t[a] + A.wrap(g[b_] - g[a])
        return np.degrees(np.stack([t[N] for N in NSTEP], 1))
    X = np.stack([np.array(NSTEP, float), np.ones(len(NSTEP))], 1)
    def slope(c, d, mask):
        s = np.linalg.lstsq(X, gap(c, d).T, rcond=None)[0][0]; return np.where(mask, s, np.nan)
    def R(label, f, unit=''):
        for nm, sel in SUB: A.reading('%s [%s]' % (label, nm), lambda c, sel=sel: f(c)[sel], unit)
    # R0 built-in null: N = 0 'DD on' and 'DD off' are the same circuit (no step, no echo): the gap's truth is 0
    a0 = alive[('off', 'D', 0)] & alive[('on', 'D', 0)]
    R('R0 null: phase, N=0 "on" - "off" (identical circuits)', lambda c: np.where(a0, np.degrees(A.wrap(ph(c, 'on', 'D', 0) - ph(c, 'off', 'D', 0))), np.nan), 'deg')
    R('R0 null: coherence, N=0 "on" - "off" (identical circuits)', lambda c: co(c, 'on', 'D', 0) - co(c, 'off', 'D', 0))
    R('R1 kick per readout, DD off (slope of phase_M - phase_D on N)', lambda c: slope(c, 'off', live['off']), 'deg')
    R('R2 kick per readout, DD on', lambda c: slope(c, 'on', live['on']), 'deg')
    R('R3 R2 - R1 (does the echo change the kick?)', lambda c: slope(c, 'on', live['both']) - slope(c, 'off', live['both']), 'deg')
    sgn = np.sign(np.nanmedian(slope((0, 1), 'off', live['both']))) or 1.0
    print('  sign of R1 (pooled median, carried): %+d; R3s = sign(R1) x (R2 - R1): BELOW = the echo shrinks the kick, ABOVE = it grows it' % sgn)
    R('R3s sign(R1) x (R2 - R1)', lambda c: sgn * (slope(c, 'on', live['both']) - slope(c, 'off', live['both'])), 'deg')
    for N in NSTEP: R('R4 coherence cost, DD off, N=%d (D - M)' % N, lambda c, N=N: co(c, 'off', 'D', N) - co(c, 'off', 'M', N))
    for N in NSTEP: R('R5 coherence cost, DD on, N=%d (D - M)' % N, lambda c, N=N: co(c, 'on', 'D', N) - co(c, 'on', 'M', N))
    for N in NSTEP: R('R6 R5 - R4, N=%d' % N, lambda c, N=N: (co(c, 'on', 'D', N) - co(c, 'on', 'M', N)) - (co(c, 'off', 'D', N) - co(c, 'off', 'M', N)))
    for N in NSTEP: R('R7 echo on the idle arm, N=%d (D on - D off)' % N, lambda c, N=N: co(c, 'on', 'D', N) - co(c, 'off', 'D', N))
A.RUNS['dd'] = dd
def tally():
    """For the null floor: how many pooled readings resolved (ABOVE/BELOW), how many CONSISTENT, how many CANNOT TELL (copies or n<5)."""
    k = {'ABOVE': 0, 'BELOW': 0, 'CONSISTENT': 0, 'CANNOT TELL': 0}
    for lab, chip, p in A.LOG: k['CANNOT TELL' if p is None else p[4]] += 1
    print('TALLY readings %d: ABOVE %d, BELOW %d, CONSISTENT %d, CANNOT TELL %d' % (len(A.LOG), k['ABOVE'], k['BELOW'], k['CONSISTENT'], k['CANNOT TELL'])); return k
def run_all(folder=HERE):
    A.run_all(folder); tally()
if __name__ == '__main__': run_all(sys.argv[1] if len(sys.argv) > 1 else HERE)
