"""Analysis exactly as PLAN_REACH.md. python analyze_reach.py [folder]. Reads job_reach_<chip>.json, raw_reach_<chip>.npz, props_reach_<chip>.json.
The ruler is deciding_runs/analyze.py, imported unchanged: per unit gap readout minus idle, median over units, 95 % bootstrap (2000 resamples),
ABOVE / BELOW / CONSISTENT / CANNOT TELL; copies compared by per-unit copy0 - copy1 (interval excludes 0 -> CANNOT TELL); signs carried (size = kick x
sign of the overall median); chips combined only by inverse variance with chi2/dof printed. A unit is one untouched qubit in one layout.
Distances are the ones written into the job at build time (coupling map before data); a check against the props couplers is printed beside."""
import os, sys, glob, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'deciding_runs'))
_argv = sys.argv; sys.argv = sys.argv[:1]; import analyze as A; sys.argv = _argv
BINS = (3, 4, 5, 6, 7, 8)
def blab(k): return 'd=%d' % k if k < 8 else 'd>=8'
def reach(J, shots):
    I = J.info; lays = I['layouts']; n = J.rec['n_qubits']; adj = A.adjacency(J.props, n) if J.props else None
    d = np.concatenate([np.array(l['d']) for l in lays]); lay_of = np.concatenate([np.full(len(l['S']), k) for k, l in enumerate(lays)])
    phys = np.concatenate([np.array(l['S']) for l in lays])
    if adj is not None:
        dp = np.concatenate([[A.bfs(adj, l['E']).get(q, 99) for q in l['S']] for l in lays]); nd = int(np.sum(dp != d))
        print('  distance check: %d of %d units sit at a different distance on the props couplers (readings use the build distances)' % (nd, len(d)))
    def arr(c, mode, st):
        return np.concatenate([J.mean_pm(c, 'fin', len(l['S']), layout=k, mode=mode, N=I['N'], setting=st) for k, l in enumerate(lays)])
    def kick(c):
        ph = lambda mode: np.arctan2(arr(c, mode, 'YY'), arr(c, mode, 'XX')); return np.degrees(A.wrap(ph('M') - ph('D'))) / I['N']
    # dead-clock guard: coherence of each unit in both arms at N = 8, both copies pooled; below 5 x shot noise -> left out of every phase reading
    coD = np.hypot(arr((0, 1), 'D', 'XX'), arr((0, 1), 'D', 'YY')); coM = np.hypot(arr((0, 1), 'M', 'XX'), arr((0, 1), 'M', 'YY'))
    thr = 5 / np.sqrt(2 * shots); alive = np.minimum(coD, coM) >= thr
    base = np.hypot(J.mean_pm((0, 1), 'fin', n, layout=-1, mode='D', N=0, setting='XX'), J.mean_pm((0, 1), 'fin', n, layout=-1, mode='D', N=0, setting='YY'))
    print('  coherence: N=0 baseline median %.3f (all %d qubits); N=8 idle median %.3f, readout median %.3f; guard threshold %.4f (5 x 1/sqrt(%d))'
          % (np.median(base), n, np.median(coD), np.median(coM), thr, 2 * shots))
    print('  dead-clock guard, units LEFT OUT per bin: ' + ', '.join('%s %d of %d' % (blab(k), np.sum(~alive & (np.minimum(d, 8) == k)), np.sum(np.minimum(d, 8) == k)) for k in BINS)
          + ', d<=2 %d of %d' % (np.sum(~alive & (d <= 2)), np.sum(d <= 2)))
    K = lambda c: np.where(alive, kick(c), np.nan)
    print('  units per bin (distinct physical qubits): ' + ', '.join('%s %d (%d)' % (blab(k), np.sum(alive & (np.minimum(d, 8) == k)), len(set(phys[alive & (np.minimum(d, 8) == k)]))) for k in BINS))
    A.reading('RR0 kick per readout, all untouched units, all layouts', K, 'deg')
    kp = K((0, 1)); sgn = np.sign(np.nanmedian(kp)) or 1.0
    print('  sign of the kick (carried, never folded): %+d; size = kick x sign.  [beside: sign of the d<=2 median %+d]' % (sgn, np.sign(np.nanmedian(kp[d <= 2])) or 0))
    sel = {k: (np.minimum(d, 8) == k) for k in BINS}; res = {}
    for k in BINS: res[k] = A.reading('RR1 kick size at %s' % blab(k), lambda c, s=sel[k]: sgn * K(c)[s], 'deg')
    v0, v1, vp = sgn * K((0,)), sgn * K((1,)), sgn * K((0, 1))
    for a, b_ in zip(BINS[:-1], BINS[1:]):
        lab = 'RR2 size(%s) - size(%s)' % (blab(a), blab(b_))
        g = A.boot2((v0 - v1)[sel[a]], (v0 - v1)[sel[b_]])
        if res[a] is None or res[b_] is None or g[3] in ('ABOVE', 'BELOW'):
            print('  %-58s copies shifted -> CANNOT TELL   (copy gap of the difference %+.4g [%+.4g, %+.4g])' % (lab, g[0], g[1], g[2])); A.LOG.append((lab, A.CUR[0], None)); continue
        m, lo, hi, o = A.boot2(vp[sel[a]], vp[sel[b_]]); p = (m, lo, hi, int(min(sel[a].sum(), sel[b_].sum())), o); A.LOG.append((lab, A.CUR[0], p))
        print('  %-58s %+.4g [%+.4g, %+.4g] deg -> %s   (copy gap %+.4g [%+.4g, %+.4g])' % (lab, m, lo, hi, o, g[0], g[1], g[2]))
    out = [(k, res[k][4] if res[k] is not None else 'CANNOT TELL', res[k]) for k in BINS]; print('  ' + rr3(out))
    print('  context, not a reading: per-layout median size per bin: ' + ' | '.join('centre %d: ' % l['centre'] + ' '.join('%s %+.2f' % (blab(k), np.nanmedian(vp[sel[k] & (lay_of == i)]) if np.any(sel[k] & (lay_of == i)) else np.nan) for k in BINS) for i, l in enumerate(lays)))
    return dict(sign=sgn, bins=out)
def rr3(out, tag='this chip'):
    """RR3, the reach EDGE: the largest bin b such that every bin from d=3 up to b is ABOVE or CANNOT TELL and b itself is ABOVE (an unbroken run of
    resolved kick from d=3; a CANNOT TELL bin neither extends nor breaks it). The first bin past the edge and its interval are printed as the bound.
    ABOVE bins beyond the first CONSISTENT/BELOW bin are listed as ISOLATED, never as the edge (a single bin out of six resolves ~2.5 % of the time
    on each side under no kick)."""
    edge, broke, iso = None, None, []
    for k, o, p in out:
        if broke is None:
            if o == 'ABOVE': edge = k
            elif o in ('CONSISTENT', 'BELOW'): broke = (k, o, p)
        elif o == 'ABOVE': iso.append(k)
    s = 'RR3 %s: reach edge (unbroken ABOVE from d=3) %s' % (tag, blab(edge) if edge else 'none')
    if broke: s += '; first bin past it %s %s, size %+.3f [%+.3f, %+.3f] deg (the bound)' % (blab(broke[0]), broke[1], broke[2][0], broke[2][1], broke[2][2])
    elif edge == 8: s += '; every bin to d>=8 resolved: the kick reaches past d=8'
    if iso: s += '; ISOLATED ABOVE beyond the edge: %s' % ', '.join(blab(k) for k in iso)
    return s + '   (bins in order: %s)' % ', '.join('%s %s' % (blab(k), o) for k, o, _ in out)
def combine():
    """As deciding_runs/analyze.summary: chips by value and width, sigma = 95 % width / 3.92, inverse-variance, chi2/dof. Returns {label: (m, lo, hi, outcome)}."""
    print('=' * 100); print('ACROSS CHIPS (values and widths; a chip whose copies shifted is listed and left out of the combination)')
    by, comb = {}, {}
    for lab, chip, p in A.LOG: by.setdefault(lab, []).append((chip, p))
    for lab, rows in by.items():
        ok = [(c, p) for c, p in rows if p is not None and np.isfinite(p[1])]
        cells = ' | '.join('%s %+.4g [%+.4g, %+.4g]' % (c.replace('ibm_', ''), p[0], p[1], p[2]) if p is not None else '%s copies shifted' % c.replace('ibm_', '') for c, p in rows)
        if len(ok) < 2: print('  %s\n      %s' % (lab, cells)); comb[lab] = None; continue
        v = np.array([p[0] for _, p in ok]); s = np.array([(p[2] - p[1]) / 3.92 for _, p in ok]); s = np.where(s > 0, s, np.nan)
        w = 1 / s ** 2; m = np.nansum(w * v) / np.nansum(w); se = 1 / np.sqrt(np.nansum(w)); lo, hi = m - 1.96 * se, m + 1.96 * se
        chi = np.nansum(((v - m) / s) ** 2) / (len(ok) - 1); o = A.outcome(lo, hi, 99); comb[lab] = (m, lo, hi, o)
        print('  %s\n      %s\n      combined %+.4g [%+.4g, %+.4g] -> %s | chip-to-chip spread chi2/dof %.2f (dof %d)' % (lab, cells, m, lo, hi, o, chi, len(ok) - 1))
    out = [(k, c[3] if c else 'CANNOT TELL', c) for k, c in ((k, comb.get('RR1 kick size at %s' % blab(k))) for k in BINS)]; print('  ' + rr3(out, 'combined'))
    return comb
def run_all(folder=HERE):
    A.LOG.clear(); per = {}
    for jf in sorted(glob.glob(os.path.join(folder, 'job_reach_*.json'))):
        chip = os.path.basename(jf)[len('job_reach_'):-5]
        if not os.path.exists(os.path.join(folder, 'raw_reach_%s.npz' % chip)): print('=' * 100); print('reach', chip, 'NOTHING RETURNED'); continue
        J = A.Job('reach', chip, folder); print('=' * 100); print('reach', chip, 'job', J.rec['job_id'], 'plan', J.rec.get('plan_sha256', '-')[:12]); A.CUR[0] = chip
        per[chip] = reach(J, J.rec['shots'])
    return per, combine()
if __name__ == '__main__': run_all(sys.argv[1] if len(sys.argv) > 1 else HERE)
