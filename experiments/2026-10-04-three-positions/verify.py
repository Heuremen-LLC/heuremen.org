#!/usr/bin/env python3
"""verify.py -- recompute every published number from the published data.

    python3 verify.py

No account, no network, no dependencies beyond the standard library. It reads
data/clocks/*.json (524 spectator-qubit readings across 7 processors, 34 KB)
and recomputes each figure, printing it beside what the page claims.

Every number it checks is on https://heuremen.org
"""
import json, math, glob, os, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
CUT = 4.0          # deg/us; spectators below this give ill-conditioned fits
PUBLISHED = {
    'pooled':      (1.0207, 0.0092),
    # 11 Oct 2026: the reading ABOVE keeps every spectator above the 4 deg/us cut, including 30
    # that report a NEGATIVE elapsed time. The reading below drops only those -- a physical bound,
    # not a tuned threshold, but chosen after the data were seen -- so the reading above is the one reported.
    'pooled_physical': (1.0331, 0.0083),
    'chi2_per_dof_physical': 0.88,
    'n_negative': 30,
    'n_above_cut': 284,
    'chi2_per_dof': 0.52,
    'p_constant':   0.79,
    'per_chip': {'ibm_fez':1.0022,'ibm_phoenix':1.0166,'ibm_aachen':1.0172,
                 'ibm_boston':1.0256,'ibm_pittsburgh':1.0266,'ibm_marrakesh':1.0282,
                 'ibm_kingston':1.0457},
}

def chi2_sf(x, k):
    if k % 2 == 0:
        t = math.exp(-x/2); s = t
        for i in range(1, k//2): t *= x/(2*i); s += t
        return min(1.0, s)
    z = math.sqrt(x); s = math.erfc(z/math.sqrt(2))
    t = math.sqrt(2/math.pi)*z*math.exp(-x/2)
    for i in range(1, (k-1)//2 + 1): s += t; t *= x/(2*i+1)
    return min(1.0, s)

def spearman(x, y):
    rx = {v:i for i,v in enumerate(sorted(x))}; ry = {v:i for i,v in enumerate(sorted(y))}
    a = [rx[v] for v in x]; b = [ry[v] for v in y]
    ma, mb = st.mean(a), st.mean(b)
    num = sum((p-ma)*(q-mb) for p,q in zip(a,b))
    den = math.sqrt(sum((p-ma)**2 for p in a) * sum((q-mb)**2 for q in b))
    return num/den if den else float('nan')

def ok(got, want, tol):
    return 'OK  ' if abs(got-want) <= tol else 'DIFF'

def main():
    files = sorted(glob.glob(os.path.join(HERE, 'data', 'clocks', '*.json')))
    if not files:
        print('no data/clocks/*.json found'); return 1
    chips = {}
    for f in files:
        d = json.load(open(f)); chips[d['chip']] = d

    print('\n  VERIFYING FROM PUBLISHED DATA -- %d processors, %d spectator readings\n'
          % (len(chips), sum(c['n'] for c in chips.values())))

    # ---- 1. the instrument note: weakly-detuned spectators report long
    print('  1. Weakly-detuned spectators report a LONGER look (published: -0.17 to -0.56, all seven)')
    for name in sorted(chips):
        c = chips[name]; r = [k['detuning_deg_per_us'] for k in c['clocks']]
        l = [k['reported_look_us'] for k in c['clocks']]
        s = spearman(r, l)
        flag = 'in range' if -0.60 <= s <= -0.10 else 'OUTSIDE'
        print('     %-16s spearman %+.3f   %s' % (name.replace('ibm_',''), s, flag))

    # ---- 2. per-chip ratio with the cut applied
    print('\n  2. Reported / stated, spectators above %.0f deg/us' % CUT)
    print('     %-14s%10s%10s%8s   %s' % ('chip','recomputed','published','n',''))
    vals = {}
    for name in sorted(chips):
        c = chips[name]
        sel = [k['reported_look_us'] for k in c['clocks'] if k['detuning_deg_per_us'] > CUT]
        if len(sel) < 8: continue
        v = st.median(sel) / c['stated_look_us']
        # bootstrap width, fixed seed so it reproduces exactly
        import random; random.seed(20261011)
        bs = [st.median(random.choices(sel, k=len(sel)))/c['stated_look_us'] for _ in range(400)]
        vals[name] = (v, st.pstdev(bs))
        p = PUBLISHED['per_chip'].get(name)
        print('     %-14s%10.4f%10.4f%8d   %s'
              % (name.replace('ibm_',''), v, p, len(sel), ok(v, p, 0.0005)))

    # ---- 3. the pooled constant
    print('\n  3. One constant across all seven')
    w = [1/s**2 for _, s in vals.values()]; S = sum(w)
    m = sum(v*x for (v,_),x in zip(vals.values(), w))/S
    sem = S**-0.5
    chi2 = sum(((v-m)/s)**2 for v, s in vals.values()); dof = len(vals)-1
    p = chi2_sf(chi2, dof)
    pm, ps = PUBLISHED['pooled']
    print('     pooled           %.4f +- %.4f   published %.4f +- %.4f   %s'
          % (m, sem, pm, ps, ok(m, pm, 0.0005)))
    print('     chi2/dof         %.2f                 published %.2f           %s'
          % (chi2/dof, PUBLISHED['chi2_per_dof'], ok(chi2/dof, PUBLISHED['chi2_per_dof'], 0.02)))
    print('     p(one constant)  %.2f                 published %.2f           %s'
          % (p, PUBLISHED['p_constant'], ok(p, PUBLISHED['p_constant'], 0.02)))

    # ---- 4. against no effect
    print('\n  4. Against no effect (ratio 1): %.1f sigma -- not established' % (abs(m-1)/sem))

    # ---- 5. what the cut is worth
    print('\n  5. Without the cut, the same arithmetic')
    allv = {}
    for name in sorted(chips):
        c = chips[name]
        sel = [k['reported_look_us'] for k in c['clocks']]
        import random; random.seed(20261011)
        v = st.median(sel)/c['stated_look_us']
        bs = [st.median(random.choices(sel,k=len(sel)))/c['stated_look_us'] for _ in range(400)]
        allv[name] = (v, st.pstdev(bs))
    w2 = [1/s**2 for _, s in allv.values()]; S2 = sum(w2)
    m2 = sum(v*x for (v,_),x in zip(allv.values(), w2))/S2
    chi2b = sum(((v-m2)/s)**2 for v, s in allv.values())
    print('     pooled %.4f   chi2/dof %.2f   p = %.4f'
          % (m2, chi2b/(len(allv)-1), chi2_sf(chi2b, len(allv)-1)))
    print('     -> one constant does NOT fit until the ill-conditioned spectators are removed.')
    print('        The cut is stated in advance and applied to every processor equally.\n')

    # ---- 6. the readings that are not elapsed times at all
    print('  6. The 30 spectators that report a NEGATIVE elapsed time')
    neg = [(n, k['detuning_deg_per_us'], k['reported_look_us']/c['stated_look_us'])
           for n, c in chips.items() for k in c['clocks']
           if k['detuning_deg_per_us'] > CUT and k['reported_look_us'] <= 0]
    tot = sum(1 for c in chips.values() for k in c['clocks'] if k['detuning_deg_per_us'] > CUT)
    print('     %d of %d above the cut  %s' % (len(neg), tot,
          ok(len(neg), PUBLISHED['n_negative'], 0) + ' / ' + ok(tot, PUBLISHED['n_above_cut'], 0)))
    print('     These come from a noisy estimator near its conditioning limit. The cut below was')
    print('     chosen AFTER the data were seen, so it is shown for transparency only. Dropping ONLY these:')
    vals2 = {}
    for name in sorted(chips):
        c = chips[name]
        sel = [k['reported_look_us'] for k in c['clocks']
               if k['detuning_deg_per_us'] > CUT and k['reported_look_us'] > 0]
        if len(sel) < 8: continue
        import random; random.seed(20261011)
        v = st.median(sel) / c['stated_look_us']
        bs = [st.median(random.choices(sel, k=len(sel)))/c['stated_look_us'] for _ in range(400)]
        vals2[name] = (v, st.pstdev(bs))
    w3 = [1/s**2 for _, s in vals2.values()]; S3 = sum(w3)
    m3 = sum(v*x for (v,_),x in zip(vals2.values(), w3))/S3; sem3 = S3**-0.5
    chi3 = sum(((v-m3)/s)**2 for v, s in vals2.values()); dof3 = len(vals2)-1
    pm3, ps3 = PUBLISHED['pooled_physical']
    print('     pooled           %.4f +- %.4f   published %.4f +- %.4f   %s'
          % (m3, sem3, pm3, ps3, ok(m3, pm3, 0.0005)))
    print('     chi2/dof         %.2f                 published %.2f           %s'
          % (chi3/dof3, PUBLISHED['chi2_per_dof_physical'], ok(chi3/dof3, PUBLISHED['chi2_per_dof_physical'], 0.02)))
    print('     vs 1 (no effect) %.1f sigma   (it is %.1f sigma with the negatives left in)'
          % (abs(m3-1)/sem3, abs(m-1)/sem))
    print('     -> post-hoc, so not reported.')
    print('        The pre-specified reading (2.2 sigma) is the one reported.\n')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
