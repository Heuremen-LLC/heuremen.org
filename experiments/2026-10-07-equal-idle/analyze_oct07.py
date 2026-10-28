"""Re-read of the 7 Oct 2026 equal-idle runs with the paper's committed estimators. Read-only on the raw data; outputs only in this folder.
Estimators: est.py (imports deciding_runs/analyze.py unchanged: wrap, boot, boot2, outcome, reading [copy gap first, then pooled reading, 95 % bootstrap over
qubits, 2000 resamples, ABOVE/BELOW/CONSISTENT/CANNOT TELL], bfs). Dead-clock guard as reach_run/analyze_reach.py. Data loading and hand mapping are the only new parts.
Usage: python3 -I analyze_oct07.py > results_oct07.txt"""
import os, sys, json, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from est import *
import sign_check
Q = '/Users/heuremen/bones-forge/parallax/out/qubit_fifth_run/data'
SHOTS = 400
LOOK6 = [1, 2, 3, 4, 6, 9]; LOOK3 = [1, 2, 4]
def rr(chip, a, b): return list(range(a, b + 1))

# family: chip, name, reads, roles file, pairs [(label, Mstem, Dstem, Ns, pool)], note
FAMS = [
 dict(chip='ibm_pittsburgh', name='FLOOR (corner, 12 looked)', runs=[1031, 1032, 1033], roles='floor_roles_ibm_pittsburgh.json', Ns=LOOK6, pairs=[('K vs D', 'K', 'D', 1)]),
 dict(chip='ibm_fez', name='FLOOR (corner, 12 looked)', runs=[1081, 1082, 1083], roles='floor_roles_ibm_fez.json', Ns=LOOK6, pairs=[('K vs D', 'K', 'D', 1)]),
 dict(chip='ibm_aachen', name='FLOOR (corner, 12 looked)', runs=[1091, 1092], roles='floor_roles_ibm_aachen.json', Ns=LOOK6, pairs=[('K vs D', 'K', 'D', 1)]),
 dict(chip='ibm_kingston', name='STATE-DARK FLOWERS (25 looked)', runs=[1191, 1192, 1193], roles='lookwide_roles_statedark_ibm_kingston.json', Ns=LOOK6,
      pairs=[('looked |0>: K vs D0', 'K', 'D0', 1), ('looked |1>: S1 vs D1', 'S1', 'D1', 0), ('looked half-way: SH vs DH', 'SH', 'DH', 0)]),
 dict(chip='ibm_kingston', name='SPACED LOOKS (25 looked)', runs=[1201, 1202, 1203], roles='lookwide_roles_spaced_ibm_kingston.json', Ns=LOOK3,
      pairs=[('looked |0>: K vs D0', 'K', 'D0', 1), ('looked |1>: S1 vs D1', 'S1', 'D1', 0), ('looked half-way: SH vs DH', 'SH', 'DH', 0)]),
 dict(chip='ibm_kingston', name='SPACED LOOKS, 3 us gaps (25 looked)', runs=[1201, 1202, 1203], roles='lookwide_roles_spaced_ibm_kingston.json', Ns=[2, 4],
      pairs=[('gapped |0>: G0 vs E0', 'G0', 'E0', 0), ('gapped |1>: G1 vs E1', 'G1', 'E1', 0), ('gapped half-way: GH vs EH', 'GH', 'EH', 0)]),
 dict(chip='ibm_kingston', name='WIRE PAINT (25 looked)', runs=[1211, 1212, 1213], roles='lookwide_roles_wirepaint_ibm_kingston.json', Ns=LOOK3,
      pairs=[('rotation 0: P0L vs P0D', 'P0L', 'P0D', 0), ('rotation 1: P1L vs P1D', 'P1L', 'P1D', 0), ('rotation 2: P2L vs P2D', 'P2L', 'P2D', 0)]),
 dict(chip='ibm_kingston', name='CARRIERS ACROSS CALIBRATIONS (25 looked)', runs=[1291, 1292, 1293], roles='lookwide_roles_carriers_ibm_kingston.json', Ns=LOOK3,
      pairs=[('index order: K vs D0', 'K', 'D0', 1), ('reversed order: KR vs D0', 'KR', 'D0', 0)]),
 dict(chip='ibm_kingston', name='CLOCK 81 (one clock, one looked bit)', runs=[1251, 1252, 1253], roles='lookwide_roles_clock81_ibm_kingston.json', Ns=LOOK3,
      pairs=[('neighbour |1>: M1L vs M1D', 'M1L', 'M1D', 0), ('neighbour |0>: M0L vs M0D', 'M0L', 'M0D', 0)]),
 dict(chip='ibm_kingston', name='CLOCK 81B (67 clocks, five looked sets)', runs=[1261, 1262, 1263], roles='lookwide_roles_clock81b_ibm_kingston.json', Ns=LOOK3,
      pairs=[('P:A vs PD', 'PA', 'PD', 0), ('P:B vs PD', 'PB', 'PD', 0), ('P:C vs PD', 'PC', 'PD', 0), ('Q:A vs QD', 'QA', 'QD', 0), ('Q:C vs QD', 'QC', 'QD', 0)]),
 dict(chip='ibm_pittsburgh', name='REACH MAP (two sets of 8 single sources)', runs=[1271, 1272, 1273], roles='lookwide_roles_reach_ibm_pittsburgh.json', Ns=LOOK3,
      pairs=[('set 0: S0L vs D', 'S0L', 'D', 1), ('set 1: S1L vs D', 'S1L', 'D', 1)]),
 dict(chip='ibm_pittsburgh', name='MEILONG / STATE-DARK FLOWERS (25 looked)', runs=[1281, 1282, 1283], roles='lookwide_roles_meilong_ibm_pittsburgh.json', Ns=LOOK3,
      pairs=[('looked |0>: K vs D0', 'K', 'D0', 1), ('looked |1>: S1 vs D1', 'S1', 'D1', 0), ('looked half-way: SH vs DH', 'SH', 'DH', 0)]),
]
BINS = [('all', lambda d: np.ones(len(d), bool)), ('d<=2', lambda d: d <= 2), ('d=3-4', lambda d: (d >= 3) & (d <= 4)), ('d>=5', lambda d: d >= 5), ('d>=3', lambda d: d >= 3)]

def load_run(chip, r):
    z = np.load('%s/fire/%s/run%d_results.npz' % (Q, chip, r)); meta = json.load(open('%s/fire/%s/run%d_meta.json' % (Q, chip, r)))
    regs = {}
    for k in z.files:
        b, g = k.split('__'); regs.setdefault(b, {})[int(g[1:])] = k
    fin = lambda b: np.asarray(z[regs[b][max(regs[b])]]).astype(complex)
    F0, F1 = fin('F0'), fin('F1'); c0 = F0.mean(0); d = F1.mean(0) - c0
    sep = np.abs(d) / np.sqrt((np.abs(F0 - c0) ** 2).mean(0) / 2)
    proj = lambda v: np.real((v - c0) * np.conj(d)) / np.abs(d) ** 2
    E = {}
    for b in regs:
        if b.endswith('|X'):
            stem = b[:-2]; v = {}
            for h in 'XY':
                p = proj(fin('%s|%s' % (stem, h)))
                v[h] = ((1 - 2 * (p > 0.5)).mean(0), (1 - 2 * p).mean(0))     # (bits variant [primary], mean-projection variant)
            E[stem] = v
    looked = {}
    for b in regs:
        if b.endswith('|X') and len(regs[b]) > 1:
            m = np.asarray(z[regs[b][1]]); looked[b[:-2]] = frozenset(np.where(np.abs(m).sum(0) > 0)[0].tolist())
    return dict(E=E, sep=sep, looked=looked, job=meta.get('job_id'), complete=meta.get('complete'), cal=meta.get('calibration_timestamp_at_fire'))

def ledger_audit(chip, r):
    L = json.load(open('%s/fire/%s/run%d_look_ledger.json' % (Q, chip, r))); dev = 0; n = 0
    for c in L['circuits']:
        if c['n_looks'] == 1 and not c['key'].startswith(('T', 'F', 'A')) and c['pre_waits_us'][0] and '|' in c['key']:
            dev = max(dev, abs(c['pre_waits_realized_us'][0] - c['pre_waits_us'][0]) / c['pre_waits_us'][0]); n += 1
    return dev, n

RESULTS = []        # (chip, family, pair label, reading label, p or None)
def run_family(F):
    chip = F['chip']; roles = json.load(open('%s/reads/%s' % (Q, F['roles']))); R = [load_run(chip, r) for r in F['runs']]
    wit = sorted(int(k) for k in roles['witness_of']); n = roles['n']
    adj = {q: set() for q in range(n)}
    for a, b in roles['edges']:
        adj[a].add(b); adj[b].add(a)
    sepok = np.all([r['sep'] > 3 for r in R], 0); nreads = len(R)
    copymap = {(0,): [0], (1,): list(range(1, nreads)), (0, 1): list(range(nreads))}
    print('=' * 110); print('%s | %s | runs %s | jobs %s' % (chip, F['name'], F['runs'], ', '.join(r['job'] for r in R)))
    dev = max(ledger_audit(chip, r)[0] for r in F['runs']); print('  design audit: dark/idle arms realised on the 16-dt grid; max |realised - stated| / stated over its plain-delay arms = %.2f %%' % (100 * dev))
    print('  copies (adaptation): the %d reads are separate jobs; copy 0 = read 1, copy 1 = %s; pooled = all. Shots 400 per read per circuit.' % (nreads, 'read 2' if nreads == 2 else 'reads 2+3 pooled'))
    clocks = np.array(wit)
    for label, Ms, Ds, pool in F['pairs']:
        Ns = F['Ns']
        keys = lambda stem: [stem + str(N) for N in Ns]
        if any(k not in R[0]['E'] for k in keys(Ms) + keys(Ds)): print('  %s: arms missing in data, skipped' % label); continue
        looked = {R[0]['looked'][k] for k in keys(Ms)}; ok_l = all(all(r['looked'][k] == R[0]['looked'][keys(Ms)[0]] for k in keys(Ms)) for r in R)
        L = sorted(R[0]['looked'][keys(Ms)[0]]); farr = set(roles['far'])
        inter = set(L) & set(wit)
        d = np.array([A.bfs(adj, L).get(q, 99) for q in wit])
        print('-' * 110); print('%s | %s' % (F['name'], label)); 
        print('  looked set read from the data (mid-register columns): %d bits, identical across N and reads: %s, subset of role file far: %s, overlaps a clock: %s' % (len(L), ok_l and len(looked) == 1, set(L) <= farr, bool(inter)))
        for variant, vname in ((0, 'bits (primary)'), (1, 'mean projection (robustness)')):
            def ex(copies, stem, N, v=variant):
                rs = copymap[copies]; X = np.mean([R[i]['E'][stem + str(N)]['X'][v] for i in rs], 0); Y = np.mean([R[i]['E'][stem + str(N)]['Y'][v] for i in rs], 0); return paper_xy_from_hands(X, Y)  # (<X>,<Y>)
            def ph(copies, stem, N): x, y = ex(copies, stem, N); return phase(x, y)[clocks]
            def co(copies, stem, N): x, y = ex(copies, stem, N); return np.hypot(x, y)[clocks]
            thr = 5 / np.sqrt(SHOTS * nreads)
            cmin = np.minimum(co((0, 1), Ms, Ns[-1]), co((0, 1), Ds, Ns[-1])); cmax = np.maximum(np.maximum(co((0, 1), Ms, Ns[0]), co((0, 1), Ds, Ns[0])), np.maximum(co((0, 1), Ms, Ns[-1]), co((0, 1), Ds, Ns[-1])))
            unprep = cmax > 1.25        # coherence cannot exceed 1: a value near sqrt(2) = both hands reading the same eigenstate, i.e. the qubit was never put on the equator in this arm (clock 81B prep Q)
            alive = sepok[clocks] & ~np.isin(clocks, L) & (cmin >= thr) & ~unprep
            kinds = {}
            def kick(copies, Ns_=Ns, icpt=False):
                tot = total_gap({N: ph(copies, Ms, N) for N in Ns_}, {N: ph(copies, Ds, N) for N in Ns_}, Ns_); return (slope_icpt if icpt else slope0)(tot, Ns_)
            def cohd(copies, N): return co(copies, Ms, N) - co(copies, Ds, N)
            def cohslope(copies): return sum(N * cohd(copies, N) for N in Ns) / np.sum(np.array(Ns, float) ** 2)
            mask = lambda f: (lambda c: np.where(alive, f(c), np.nan))
            N4 = 4 if 4 in Ns else Ns[-1]
            if variant == 0:
                print('  clocks: %d in role file, %d kept (cloud separation > 3 in every read, not looked at, coherence >= 5/sqrt(%d) = %.4f in both arms at N=%d [dead-clock guard as reach_run]; %d dropped as never put on the equator in this arm, coherence > 1.25)' % (len(clocks), alive.sum(), SHOTS * nreads, thr, Ns[-1], int((unprep & sepok[clocks]).sum())))
                print('  clocks per distance bin (couplings to nearest looked bit): ' + ', '.join('%s %d' % (bn, (alive & bf(d)).sum()) for bn, bf in BINS))
                A.CUR[0] = chip
                slopes = {}
                for bn, bf in BINS:
                    sel = bf(d); p = A.reading('%s kick/readout, slope vs N (N=%s) [%s]' % (label[:22], '/'.join(map(str, Ns)), bn), lambda c, sel=sel: mask(kick)(c)[sel], 'deg'); RESULTS.append((chip, F['name'], label, 'kick ' + bn, p, pool))
                rs = A.reading('%s kick, N=%d only (R1 form) [all]' % (label[:22], N4), mask(lambda c: slope0(total_gap({N4: ph(c, Ms, N4)}, {N4: ph(c, Ds, N4)}, [N4]), [N4])), 'deg'); RESULTS.append((chip, F['name'], label, 'kickN4', rs, 0))
                rs = A.reading('%s kick, slope with intercept (R2.3 form) [all]' % label[:22], mask(lambda c: kick(c, icpt=True)), 'deg'); RESULTS.append((chip, F['name'], label, 'kickicpt', rs, 0))
                for bn, bf in (BINS[0], BINS[1], BINS[4]):
                    sel = bf(d); p = A.reading('%s coherence M-D, slope vs N [%s]' % (label[:22], bn), lambda c, sel=sel: mask(cohslope)(c)[sel]); RESULTS.append((chip, F['name'], label, 'coh ' + bn, p, 0))
                p = A.reading('%s coherence M-D at N=%d [all]' % (label[:22], N4), mask(lambda c: cohd(c, N4))); RESULTS.append((chip, F['name'], label, 'cohN4', p, 0))
                # per-read medians (context, not readings)
                def ph_read(i, stem, N):
                    x = paper_xy_from_hands(R[i]['E'][stem + str(N)]['X'][0], R[i]['E'][stem + str(N)]['Y'][0]); return phase(*x)[clocks]
                pr = []
                for i in range(nreads):
                    kk = slope0(total_gap({N: ph_read(i, Ms, N) for N in Ns}, {N: ph_read(i, Ds, N) for N in Ns}, Ns), Ns); pr.append(np.nanmedian(np.where(alive, kk, np.nan)))
                print('  context, per-read median kick (all kept clocks, not a reading): ' + ' | '.join('read %d %+.3f' % (i + 1, v) for i, v in enumerate(pr)))
                if alive.sum() < 5:
                    for i in range(nreads):
                        kk = slope0(total_gap({N: ph_read(i, Ms, N) for N in Ns}, {N: ph_read(i, Ds, N) for N in Ns}, Ns), Ns)
                        print('  single-clock values, read %d: clock(s) %s kick %s deg/readout, distance %s' % (i + 1, clocks[alive].tolist(), np.round(kk[alive], 2).tolist(), d[alive].tolist()))
                POOL.setdefault(chip, []).append(dict(label=F['name'] + ' / ' + label, kick=mask(kick), d=d, pool=pool, bins=None))
            else:
                p = A.reading('%s kick, mean-projection variant [all]' % label[:22], mask(kick), 'deg'); RESULTS.append((chip, F['name'], label, 'kick_meanproj', p, 0))
                p2 = A.reading('%s kick, mean-projection variant [d<=2]' % label[:22], lambda c: mask(kick)(c)[d <= 2], 'deg')
                p3 = A.reading('%s kick, mean-projection variant [d>=3]' % label[:22], lambda c: mask(kick)(c)[d >= 3], 'deg')
POOL = {}
def chip_pools():
    print('=' * 110); print('PER-CHIP POOL (units = clock x family; plain |0> back-to-back arms only: the K-type arms with the looked bits in |0>; clocks repeat across families on a chip, so the pooled n counts repeated measures)')
    print('  Used here only to give one per-chip line by distance. Per-family readings above are the readings of record.')
    for chip, lst in POOL.items():
        use = [u for u in lst if u['pool']]
        if not use: continue
        print('-' * 110); print('%s pool of: %s' % (chip, '; '.join(u['label'] for u in use)))
        for bn, bf in BINS:
            f = lambda c, bn=bn, bf=bf: np.concatenate([u['kick'](c)[bf(u['d'])] for u in use])
            p = A.reading('%s POOL kick/readout [%s]' % (chip.replace('ibm_', ''), bn), f, 'deg'); RESULTS.append((chip, 'POOL', 'plain |0> arms', 'kick ' + bn, p, 0))

if __name__ == '__main__':
    print('RE-READ OF THE 7 OCT 2026 EQUAL-IDLE RUNS WITH THE PAPER\'S COMMITTED ESTIMATORS (read-only; no IBM jobs)')
    print('Sign: paper convention, phase = atan2(<Y>,<X>); kick = (readout - idle phase)/N, so the 21-22 Oct kick is NEGATIVE; BELOW 0 = same direction as 21-22 Oct.'); print()
    ok = sign_check.run(); print(); 
    for F in FAMS: run_family(F)
    chip_pools()
    json.dump([(c, f, l, k, p if p is None else [float(x) if not isinstance(x, str) else x for x in p]) for c, f, l, k, p, _ in RESULTS], open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_oct07.json'), 'w'))
