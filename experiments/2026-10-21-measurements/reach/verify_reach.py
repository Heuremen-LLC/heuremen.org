#!/usr/bin/env python3
"""verify_reach.py -- recompute the "reach" quantities from the per-spectator data.

    python3 verify_reach.py

No account, no network, no dependencies beyond the standard library. It reads data/*.json
(layout_A_pool.json, layout_B.json, clusters.json, floor.json; about 180 KB in all) and recomputes
how far the effect of a readout reaches:

  1. fall-off with distance (3-5 / 6-8 / 9+ spacings, idle control)
  2. three predictors ranked against each spectator's reading
  3. dose-response by number of measured qubits within two spacings (layout A, layout B) and the
     change per additional measured qubit
  4. four clusters of readouts around four central spectators

Nothing here is tuned: every rule (which spectators, which estimator, how bins are formed) is named
where it is applied.
"""
import json, math, os, random, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
D = lambda f: json.load(open(os.path.join(HERE, 'data', f)))
NAN = float('nan')
fin = lambda x: x is not None and not math.isnan(x)

# column positions in the per-spectator rows of layout_A_pool / layout_B / clusters
C_ID, C_N2, C_INT, C_NEAR, C_RATE, C_FIT, C_39 = range(7)
# column positions in the rows of floor.json
F_ID, F_DIST, F_D, F_K, F_K2 = range(5)


def med(v):
    v = [x for x in v if fin(x)]
    return st.median(v) if v else NAN


def mean(v):
    v = [x for x in v if fin(x)]
    return sum(v) / len(v) if v else NAN


def ranks(x):
    """average ranks (ties share the mean rank)."""
    o = sorted(range(len(x)), key=lambda i: x[i]); r = [0.0] * len(x); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and x[o[j + 1]] == x[o[i]]: j += 1
        for k in range(i, j + 1): r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(a, b):
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    num = sum((p - ma) * (q - mb) for p, q in zip(a, b))
    den = math.sqrt(sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b))
    return num / den if den else NAN


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def ols_slope(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)


# =====================================================================================================
def section_fall_off():
    print('  1. Fall-off with distance.  Ratio of reported/stated (measured arm) to reported/laid (idle arm), by spacing band.')
    F = D('floor.json'); bands = [('3-5', 3.0, 5.5), ('6-8', 5.5, 8.5), ('9+', 8.5, 1e9)]

    def band_table(chip, reads_idx):
        c = [x for x in F['chips'] if x['chip'] == chip][0]; rd = [c['reads'][i] for i in reads_idx if i < len(c['reads'])]
        out = {}
        for name, lo, hi in bands + [('all', 3.0, 1e9)]:
            per = {'D': [], 'K': [], 'K2': []}
            for r in rd:
                rows = [x for x in r['clocks'] if lo <= x[F_DIST] < hi]
                per['D'].append(med([x[F_D] for x in rows]) / c['idle_us_per_readout'])
                per['K'].append(med([x[F_K] for x in rows]) / c['stated_measure_us'])
                per['K2'].append(med([x[F_K2] for x in rows]) / c['stated_measure_2_us'])
            out[name] = {k: mean(v) for k, v in per.items()}; out[name]['n'] = len([x for x in rd[0]['clocks'] if lo <= x[F_DIST] < hi]); out[name]['reads'] = per
        return c, out

    c, t12 = band_table('ibm_pittsburgh', [0, 1]); _, t123 = band_table('ibm_pittsburgh', [0, 1, 2]); _, t1 = band_table('ibm_pittsburgh', [0])
    runs = [r['run'] for r in c['reads']]
    print('     data: ibm_pittsburgh, 7 Oct 2026, runs %s; 73 spectators 3 or more spacings from every measured qubit.' % runs)
    print('     Three arms in one job per read: idle (identical idle time, nothing measured), measure (the chip\'s own %.2f us readout), measure_2 (a %.2f us readout).'
          % (c['stated_measure_us'], c['stated_measure_2_us']))
    print('     ratio = (reported/stated on the measured arm) / (reported/laid on the idle arm); reads 1 and 2 of the job in the first table.')
    print('\n     %-34s %8s' % ('measure / idle, reads 1+2', 'value'))
    for name in ('3-5', '6-8', '9+'):
        v = t12[name]['K'] / t12[name]['D']; print('       %-32s %8.3f' % ('%s spacings (n=%d)' % (name, t12[name]['n']), v))
    v = t12['all']['D']; print('       %-32s %8.3f   <- idle arm, absolute (reported/laid)' % ('idle control, all 73', v))
    print('\n     The same quantity as more of the data is admitted (measure / idle):')
    print('       %-34s %8s %8s %8s' % ('', '3-5', '6-8', '9+'))
    rows = [('Pittsburgh, read 1 only', t1), ('Pittsburgh, reads 1+2', t12), ('Pittsburgh, all 3 reads', t123)]
    for lab, t in rows: print('       %-34s %8.3f %8.3f %8.3f' % (lab, *(t[b]['K'] / t[b]['D'] for b in ('3-5', '6-8', '9+'))))
    print('       %-34s %8.3f %8.3f %8.3f   (measure_2 arm, all 3 reads)' % ('', *(t123[b]['K2'] / t123[b]['D'] for b in ('3-5', '6-8', '9+'))))
    print('       %-34s %8s %8s %8s' % ('read-to-read, measure/idle:', '', '', ''))
    for i in range(3):
        print('         read %d (run %d)                  %8.3f %8.3f %8.3f' % (i + 1, runs[i], *(t123[b]['reads']['K'][i] / t123[b]['reads']['D'][i] for b in ('3-5', '6-8', '9+'))))
    for chip in ('ibm_fez', 'ibm_aachen'):
        cc, tt = band_table(chip, [0, 1, 2]); cc1, t_1 = band_table(chip, [0])
        print('       %-34s %8.3f %8.3f %8.3f   (runs %s; read 1 alone: %.3f %.3f %.3f)'
              % (chip + ', all reads', *(tt[b]['K'] / tt[b]['D'] for b in ('3-5', '6-8', '9+')), [r['run'] for r in cc['reads']], *(t_1[b]['K'] / t_1[b]['D'] for b in ('3-5', '6-8', '9+'))))
    # is the near band different from the far band? (spread of the three reads, in band-ratio units)
    sd = st.pstdev([t123['3-5']['reads']['K'][i] / t123['3-5']['reads']['D'][i] for i in range(3)])
    print('     -> Over all three reads of the Pittsburgh job the near-band excess is %+.3f and the 6-8 band is'
          % (t123['3-5']['K'] / t123['3-5']['D'] - 1))
    print('        larger than the 3-5 band; the read-to-read spread of the 3-5 band is %.3f, as large as the effect. The two other chips run the same design show no near-band excess.' % sd)
    print('     NOTE the idle-control value is an absolute reported/laid value (idle arm), while the band ratios already divide by it; the two are not on one scale.')
    return t12, t123


# =====================================================================================================
def pool_rows(pool, skip_phoenix=False):
    """one dict per (spectator, read) in layout A: reading divided by the stated readout length, 0.3 < ratio < 5 kept."""
    out = []
    for cd in pool['chip_days']:
        if skip_phoenix and cd['chip'] == 'ibm_phoenix': continue
        for r in cd['reads']:
            for x in r['clocks']:
                y = x[C_FIT] / cd['stated_readout_us'] if fin(x[C_FIT]) else NAN
                out.append(dict(chip=cd['chip'], day=cd['date'], run=r['run'], clock=x[C_ID], n2=x[C_N2], inten=x[C_INT], near=x[C_NEAR], rate=x[C_RATE], y=y))
    return out


def boot_rho(rows, key, nb=1000, seed=20261021):
    rng = random.Random(seed); n = len(rows); y = [r['y'] for r in rows]; xs = {k: [r[k] for r in rows] for k in ('n2', 'inten', 'near')}
    d1, d2 = [], []
    for _ in range(nb):
        idx = [rng.randrange(n) for _ in range(n)]
        yy = [y[i] for i in idx]; rc = spearman([xs['n2'][i] for i in idx], yy); ri = spearman([xs['inten'][i] for i in idx], yy); rn = spearman([xs['near'][i] for i in idx], yy)
        d1.append(ri - rc); d2.append(ri - abs(rn))
    q = lambda v: (sorted(v)[int(0.025 * nb)], sorted(v)[int(0.975 * nb)])
    return q(d1), q(d2)


def section_predictors(pool, clus):
    print('\n  2. Three predictors ranked against each spectator\'s reading (Spearman rank correlation with reported/stated).')
    print('     Predictors: intensity (sum of 1/r^2 over all measured qubits), count of measured qubits within two spacings, nearest measured-qubit distance.')
    sets = []
    sets.append(('layout A pool (5 chips, 8 chip-days, no Phoenix)', [r for r in pool_rows(pool, True) if fin(r['y']) and 0.3 < r['y'] < 5]))
    fez = [cd for cd in clus['chip_days'] if cd['chip'] == 'ibm_fez'][0]

    def fez_rows(reads):
        out = []
        for r in fez['reads']:
            if r['run'] not in reads: continue
            for x in r['clocks']:
                y = x[C_FIT] / fez['stated_readout_us'] if fin(x[C_FIT]) else NAN
                if fin(y) and 0.3 < y < 5: out.append(dict(n2=x[C_N2], inten=x[C_INT], near=x[C_NEAR], y=y))
        return out
    sets.append(('four-cluster chip (Fez), read 1 only', fez_rows({921})))
    sets.append(('four-cluster chip (Fez), all 3 reads', fez_rows({921, 922, 923})))
    print('\n     %-50s %5s | %8s %8s %8s' % ('', 'n', 'count', 'intens.', 'nearest'))
    for lab, rows in sets:
        y = [r['y'] for r in rows]
        rc, ri, rn = (spearman([r[k] for r in rows], y) for k in ('n2', 'inten', 'near'))
        print('     %-50s %5d | %+8.3f %+8.3f %+8.3f' % (lab, len(rows), rc, ri, rn))
        best = max((abs(rc), 'count'), (abs(ri), 'intensity'), (abs(rn), 'nearest'))[1]
        (lo, hi), (lo2, hi2) = boot_rho(rows, None)
        print('        best by |rho|: %s.  bootstrap 95%% interval (resampling spectator-reads, 1000 draws) of rho_intensity - rho_count: [%+.3f, %+.3f];  of rho_intensity - |rho_nearest|: [%+.3f, %+.3f]'
              % (best, lo, hi, lo2, hi2))
    print('     -> intensity has the largest |rho| in each set. In the pool all three correlations are weak (|rho| <= %.2f, i.e. under 3%% of rank variance)'
          % max(abs(spearman([r[k] for r in sets[0][1]], [r['y'] for r in sets[0][1]])) for k in ('n2', 'inten', 'near')))
    print('        and the intervals above show whether the margin over the other two predictors is distinguishable from zero.')
    return sets


# =====================================================================================================
def chipday_bins(cd, minrate=None):
    """per read: median over spectators in each bin (0,1,2,>=3 measured qubits within two spacings; Phoenix: within two links) of reported/stated; mean over reads."""
    per = [[] for _ in range(4)]
    for r in cd['reads']:
        for b in range(4):
            v = [x[C_FIT] for x in r['clocks'] if min(x[C_N2], 3) == b and fin(x[C_FIT]) and (minrate is None or abs(x[C_RATE]) > minrate)]
            all_v = [x[C_FIT] for x in r['clocks'] if min(x[C_N2], 3) == b and (minrate is None or abs(x[C_RATE]) > minrate)]
            if all_v: per[b].append(med([x[C_FIT] for x in r['clocks'] if min(x[C_N2], 3) == b and (minrate is None or abs(x[C_RATE]) > minrate)]) / cd['stated_readout_us'])
    return [mean(p) for p in per], per


def section_dose(pool, B):
    print('\n  3. Dose-response: spectators binned by number of measured qubits within two spacings (0, 1, 2, >=3).')
    chips = {}
    for cd in pool['chip_days']:
        b, per = chipday_bins(cd); chips.setdefault(cd['chip'], []).append((cd, b))
    perchip = {c: [mean([b[i] for _, b in v]) for i in range(4)] for c, v in chips.items()}
    A = [mean([perchip[c][i] for c in perchip]) for i in range(4)]
    print('\n     Layout A pools %d processors (%d chip-days, %d reads): the mean over processors of each processor\'s mean of per-read medians.'
          % (len(perchip), sum(len(v) for v in chips.values()), sum(len(cd['reads']) for v in chips.values() for cd, _ in v)))
    print('     %-34s %7s %7s %7s %7s' % ('', '0', '1', '2', '>=3'))
    print('     %-34s %7.3f %7.3f %7.3f %7.3f' % ('layout A (mean over 6 processors)', *A))
    for c in sorted(perchip): print('       %-32s %7.3f %7.3f %7.3f %7.3f   (%s)' % (c.replace('ibm_', ''), *perchip[c], ', '.join(cd['date'][5:] + ' x%d' % len(cd['reads']) for cd, _ in chips[c])))
    up = sum(1 for c in perchip if perchip[c][3] > perchip[c][0])
    print('       processors where >=3 exceeds 0: %d of %d.  Monotone 0<1<2<3+ on: %s'
          % (up, len(perchip), ', '.join(c.replace('ibm_', '') for c in sorted(perchip) if perchip[c][0] < perchip[c][1] < perchip[c][2] < perchip[c][3]) or 'none'))
    # sensitivity: drop spectators precessing under 4 deg/us (the cut the duration table uses)
    chips4 = {}
    for cd in pool['chip_days']:
        b, _ = chipday_bins(cd, 4.0); chips4.setdefault(cd['chip'], []).append(b)
    pc4 = {c: [mean([b[i] for b in v]) for i in range(4)] for c, v in chips4.items()}
    A4 = [mean([pc4[c][i] for c in pc4]) for i in range(4)]
    n4 = [sum(1 for cd in pool['chip_days'] for r in cd['reads'] for x in r['clocks'] if min(x[C_N2], 3) == b and fin(x[C_FIT]) and abs(x[C_RATE]) > 4.0) for b in range(4)]
    n0 = [sum(1 for cd in pool['chip_days'] for r in cd['reads'] for x in r['clocks'] if min(x[C_N2], 3) == b and fin(x[C_FIT])) for b in range(4)]
    print('     %-34s %7.3f %7.3f %7.3f %7.3f   (only spectators with |rate| > 4 deg/us, the cut the duration table uses; spectator-reads kept %s of %s)' % ('same, rate cut applied', *A4, '/'.join(map(str, n4)), '/'.join(map(str, n0))))

    # ---- layout B
    print('\n     Layout B (an independent second random draw of which qubits are measured), 7 Oct:')
    byc = {cd['chip']: cd for cd in B['chip_days']}
    pitt = byc['ibm_pittsburgh']; bp, per = chipday_bins(pitt)
    print('       ibm_pittsburgh layout B runs %s: %s' % ([r['run'] for r in pitt['reads']], ' '.join('%.3f' % x for x in bp)))
    print('         0 measured in range:  %.3f' % bp[0])
    print('         >=3 measured in range: %.3f' % bp[3])
    pa = [(cd, b) for cd, b in chips['ibm_pittsburgh']]
    print('       ibm_pittsburgh is also one of the six processors pooled in layout A (layout A there, 4 Oct and 5 Oct: %s).' % '; '.join(' '.join('%.3f' % x for x in b) for _, b in pa))
    print('         On that same processor, layout A on 4 Oct gives >=3 -> %.3f (0 -> %.3f); layout B three days later gives %.3f -> %.3f.' % (pa[0][1][3], pa[0][1][0], bp[0], bp[3]))
    mar = byc['ibm_marrakesh']; bm, _ = chipday_bins(mar)
    ma = chips['ibm_marrakesh'][0][1]
    print('       ibm_marrakesh layout B runs %s: %s   (layout A, same processor, 4 Oct: %s)' % ([r['run'] for r in mar['reads']], ' '.join('%.3f' % x for x in bm), ' '.join('%.3f' % x for x in ma)))
    print('         -> layout B on ibm_marrakesh gives >=3 -> %.3f, below its zero bin %.3f.' % (bm[3], bm[0]))
    print('         Marrakesh layout B has %d spectators with >=3 measured qubits in range; read-to-read: %s' % (sum(1 for x in mar['reads'][0]['clocks'] if x[C_N2] >= 3), ' / '.join('%.3f' % (med([x[C_FIT] for x in r['clocks'] if x[C_N2] >= 3]) / mar['stated_readout_us']) for r in mar['reads'])))

    # ---- change per additional measured qubit
    print('\n     Change in reading per additional measured qubit in range:')
    print('       table, (>=3 bin - 0 bin) / 3 rungs (layout A)                      %+.3f' % ((A[3] - A[0]) / 3))
    print('       table, layout B Pittsburgh, (>=3 - 0) / 3                          %+.3f' % ((bp[3] - bp[0]) / 3))
    rows = [(x[C_N2], x[C_FIT] / cd['stated_readout_us']) for cd in pool['chip_days'] if cd['chip'] != 'ibm_phoenix' for r in cd['reads'] for x in r['clocks'] if fin(x[C_FIT]) and 0.3 < x[C_FIT] / cd['stated_readout_us'] < 5]
    print('       least-squares slope, every spectator-read in layout A, 0..4 in range  %+.3f   (n=%d)' % (ols_slope([a for a, _ in rows], [b for _, b in rows]), len(rows)))
    rowsb = [(x[C_N2], x[C_FIT] / pitt['stated_readout_us']) for r in pitt['reads'] for x in r['clocks'] if fin(x[C_FIT]) and 0.3 < x[C_FIT] / pitt['stated_readout_us'] < 5]
    print('       least-squares slope, every spectator-read in layout B Pittsburgh      %+.3f   (n=%d)' % (ols_slope([a for a, _ in rowsb], [b for _, b in rowsb]), len(rowsb)))
    # same-clock change between layout A (4 Oct) and layout B (7 Oct) on one processor
    def shared_change(cdA, cdB):
        def tab(cd):
            t = []
            for r in cd['reads']:
                t.append({x[C_ID]: (x[C_FIT] / cd['stated_readout_us'] if fin(x[C_FIT]) else NAN, x[C_N2]) for x in r['clocks']})
            return t
        TA, TB = tab(cdA), tab(cdB)
        sh = [w for w in TA[0] if all(w in t for t in TA + TB)]
        mn = lambda v: sum(v) / len(v)                                   # NaN propagates
        dR = {w: mn([t[w][0] for t in TB]) - mn([t[w][0] for t in TA]) for w in sh}; ds = {w: TB[-1][w][1] - TA[-1][w][1] for w in sh}
        return sh, dR, ds
    for chip in ('ibm_pittsburgh', 'ibm_marrakesh'):
        sh, dR, ds = shared_change(chips[chip][0][0], byc[chip])
        k1 = [w for w in sh if ds[w] == 1]; okw = [w for w in sh if fin(dR[w])]
        print('       %s: %d spectators present in both layouts: least-squares change in reading per extra measured qubit in range %+.3f;  median change of the %d that gained exactly one: %+.3f'
              % (chip.replace('ibm_', ''), len(sh), ols_slope([ds[w] for w in okw], [dR[w] for w in okw]), len(k1), med([dR[w] for w in k1])))
    sh, dR, ds = shared_change(chips['ibm_pittsburgh'][0][0], pitt)
    plus1 = med([dR[w] for w in sh if ds[w] == 1])
    print('       Pittsburgh median change for spectators that gained exactly one measured qubit: %+.3f (%d spectators); this is a subgroup median, not a slope.' % (plus1, sum(1 for w in sh if ds[w] == 1)))
    return A, bp, bm


# =====================================================================================================
def section_clusters(clus):
    print('\n  4. Four clusters of measured qubits around four central spectators.')
    c = [x for x in clus['chip_days'] if x['chip'] == 'ibm_fez'][0]; O = c['stated_readout_us']
    runs = [r['run'] for r in c['reads']]
    print('     data: ibm_fez, 7 Oct 2026, runs %s (same layout, three consecutive jobs); stated readout %.3f us; centres (spectator qubit numbers) %s; %d measured qubits.'
          % (runs, O, c['centres'], c['n_measured_qubit_bits']))
    print('\n     reported/stated at each centre (nan = the spectator\'s phase signal was too weak to fit):')
    print('       %-10s %s' % ('centre', '  '.join('read %d (%d)' % (i + 1, r) for i, r in enumerate(runs))))
    vals = {}
    for w in c['centres']:
        row = []
        for r in c['reads']:
            x = [x for x in r['clocks'] if x[C_ID] == w]
            row.append(x[0][C_FIT] / O if x and fin(x[0][C_FIT]) else NAN)
        vals[w] = row
        print('       %-10s %s' % (w, '  '.join('%10.3f' % v if fin(v) else '       nan' for v in row)))
    allv = [v for row in vals.values() for v in row if fin(v)]
    print('       readable: %d of %d (centre, read) cells; values %s; median of readable %.3f' % (len(allv), 4 * len(runs), ', '.join('%.3f' % v for v in sorted(allv)), med(allv)))
    print('       The 4th centre (131) was unreadable on reads 1-2 and read %.3f on read 3.' % (vals[131][2] if 131 in vals else NAN))
    print('       read 3 of the same job: %s.  Per-read mean of the readable centres: %s'
          % (', '.join('%d: %.3f' % (w, vals[w][2]) for w in vals if fin(vals[w][2])), ' / '.join('%.3f' % mean([vals[w][i] for w in vals]) for i in range(len(runs)))))
    print('       %d of the readable centre-values are below 1 or negative: %s (a negative elapsed time is not a duration; this is the estimator failing near its limit).'
          % (sum(1 for v in allv if v < 1.0), ', '.join('%.3f' % v for v in allv if v < 1.0)))
    # the open field
    print('\n     Open field: the same job\'s spectators with no measured qubit within two spacings (median of reported/stated):')
    for i, r in enumerate(c['reads']):
        z = [x[C_FIT] / O for x in r['clocks'] if x[C_N2] == 0 and fin(x[C_FIT])]; z39 = [x[C_39] / O for x in r['clocks'] if x[C_N2] == 0 and fin(x[C_39])]
        print('       read %d (run %d): fit estimator %.3f (n=%d) | 3-and-9 estimator %.3f' % (i + 1, r['run'], med(z), len(z), med(z39)))
    zall = [x[C_FIT] / O for r in c['reads'] for x in r['clocks'] if x[C_N2] == 0 and fin(x[C_FIT])]
    print('       three reads pooled (fit estimator): %.3f' % med(zall))
    # the same layout on the other three processors
    print('\n     The same four-cluster layout, same day, three other processors (reported/stated at the centres, reads 1/2/3):')
    for cc in clus['chip_days']:
        if cc['chip'] == 'ibm_fez': continue
        O2 = cc['stated_readout_us']; out = []
        for w in cc['centres']:
            row = []
            for r in cc['reads']:
                x = [x for x in r['clocks'] if x[C_ID] == w]; row.append(x[0][C_FIT] / O2 if x and fin(x[0][C_FIT]) else NAN)
            out.append('%d: %s' % (w, '/'.join('%.2f' % v if fin(v) else 'nan' for v in row)))
        print('       %-16s runs %s  %s' % (cc['chip'], [r['run'] for r in cc['reads']], ' | '.join(out)))
    return vals


def main():
    pool, B, clus = D('layout_A_pool.json'), D('layout_B.json'), D('clusters.json')
    print('\n  REACH OF A READOUT -- %d chip-days of layout A, 2 of layout B, %d four-cluster chips, 3 floor chips\n'
          % (len(pool['chip_days']), len(clus['chip_days'])))
    section_fall_off()
    section_predictors(pool, clus)
    section_dose(pool, B)
    section_clusters(clus)
    print()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
