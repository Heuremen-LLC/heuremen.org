"""Analysis exactly as PLAN.md. python analyze.py [folder]. Reads job_<run>_<chip>.json, raw_<run>_<chip>.npz, props_<run>_<chip>.json.
Every reading is a GAP with a 95 % bootstrap interval (2000 resamples over qubits) and one of four outcomes: ABOVE / BELOW / CONSISTENT / CANNOT TELL.
Copies: per qubit copy0 - copy1; interval contains 0 -> the pooled value is the reading; interval excludes 0 -> CANNOT TELL (copies shifted), both printed.
No thresholds in sigma; no outcome is decided by comparing two labels."""
import json, glob, os, sys
from collections import deque
import numpy as np
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)); B = 2000
RNG = np.random.default_rng(20261032)
def outcome(lo, hi, n):
    if n < 5 or not np.isfinite(lo): return 'CANNOT TELL'
    return 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def tmean(x, axis=None):
    x = np.sort(x, axis=-1); k = int(0.2 * x.shape[-1]); return np.mean(x[..., k:x.shape[-1] - k], axis=-1)
def boot(x, stat='median'):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 5: return np.nan, np.nan, np.nan, len(x), 'CANNOT TELL'
    f = (lambda y: np.median(y, -1)) if stat == 'median' else tmean
    b = f(RNG.choice(x, (B, len(x)))); lo, hi = np.percentile(b, [2.5, 97.5]); return float(f(x[None])[0]), lo, hi, len(x), outcome(lo, hi, len(x))
def boot2(a, b_):
    a = np.asarray(a, float); b_ = np.asarray(b_, float); a = a[np.isfinite(a)]; b_ = b_[np.isfinite(b_)]
    if min(len(a), len(b_)) < 5: return np.nan, np.nan, np.nan, 'CANNOT TELL'
    d = np.median(RNG.choice(a, (B, len(a))), 1) - np.median(RNG.choice(b_, (B, len(b_))), 1); lo, hi = np.percentile(d, [2.5, 97.5])
    return float(np.median(a) - np.median(b_)), lo, hi, outcome(lo, hi, min(len(a), len(b_)))
LOG = []; CUR = ['']
def reading(label, f, unit='', stat='median'):
    """f(copies) -> per-qubit array for copies (0,), (1,), (0, 1). Copy gap first, then the pooled reading."""
    v0, v1, vp = f((0,)), f((1,)), f((0, 1))
    g = boot(v0 - v1, stat)
    if g[4] == 'ABOVE' or g[4] == 'BELOW':
        print('  %-58s copies shifted %+.4g [%+.4g, %+.4g] | copy0 %+.4g copy1 %+.4g -> CANNOT TELL' % (label, g[0], g[1], g[2], np.nanmedian(v0), np.nanmedian(v1)))
        LOG.append((label, CUR[0], None)); return None
    p = boot(vp, stat); LOG.append((label, CUR[0], p))
    print('  %-58s %+.4g [%+.4g, %+.4g] %s n=%d -> %s   (copy gap %+.4g [%+.4g, %+.4g])' % (label, p[0], p[1], p[2], unit, p[3], p[4], g[0], g[1], g[2])); return p
def bits(arr, nbits):
    return np.unpackbits(np.asarray(arr, np.uint8), axis=1, bitorder='big')[:, ::-1][:, :nbits].astype(np.int8)
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
def adjacency(props, n):
    adj = {q: set() for q in range(n)}
    for g in props['gates']:
        if g['gate'] in ('cz', 'ecr', 'cx') and len(g['qubits']) == 2: a, c = g['qubits']; adj[a].add(c); adj[c].add(a)
    return adj
def bfs(adj, src):
    d = {s: 0 for s in src}; q = deque(src)
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in d: d[v] = d[u] + 1; q.append(v)
    return d
class Job:
    def __init__(self, run, chip, folder=HERE, synth=None):
        self.rec = json.load(open(os.path.join(folder, 'job_%s_%s.json' % (run, chip)))); self.meta = self.rec['meta']; self.info = self.rec['info']
        self.z = synth if synth is not None else np.load(os.path.join(folder, 'raw_%s_%s.npz' % (run, chip)))
        self.props = json.load(open(os.path.join(folder, 'props_%s_%s.json' % (run, chip))))
    def idx(self, copies, **kw): return [i for i, m in enumerate(self.meta) if m['copy'] in copies and all(m.get(k) == v for k, v in kw.items())]
    def get(self, i, reg, nbits):
        a = self.z['%s_%d' % (reg, i)]; return a if a.dtype == np.int8 else bits(a, nbits)
    def mean_pm(self, copies, reg, nbits, **kw):   # per bit <1 - 2b>, pooled over matching pubs
        x = np.concatenate([self.get(i, reg, nbits) for i in self.idx(copies, **kw)]); return np.mean(1 - 2 * x, 0)
    def p1(self, copies, reg, nbits, **kw):
        x = np.concatenate([self.get(i, reg, nbits) for i in self.idx(copies, **kw)]); return np.mean(x, 0)
# ------------------------------------------------------------------------------------------------ run 1: distance
def distance(J):
    I = J.info; n = len(I['E_ball']) + len(I['S_ball']); adj = adjacency(J.props, n)
    out = {}
    for lay in ('ball', 'spread'):
        E, S = I['E_%s' % lay], I['S_%s' % lay]; nS = len(S); dd = bfs(adj, E); d = np.array([dd.get(s, 99) for s in S])
        inv = np.array([sum(1.0 / bfs(adj, [s]).get(e, 99) ** 2 for e in E) for s in S])
        def kick(c, mo, lay=lay, nS=nS):
            ph = lambda mode: np.arctan2(J.mean_pm(c, 'fin', nS, layout=lay, mode=mode, N=8, setting='YY'), J.mean_pm(c, 'fin', nS, layout=lay, mode=mode, N=8, setting='XX'))
            return np.degrees(wrap(ph('M') - ph('D'))) / 8
        out[lay] = (d, inv, kick)
    print(' layout ball: spectators at distance <=2 / 3-4 / 5+: %s' % [int(x) for x in (np.sum(out['ball'][0] <= 2), np.sum((out['ball'][0] >= 3) & (out['ball'][0] <= 4)), np.sum(out['ball'][0] >= 5))])
    d, inv, kick = out['ball']
    k = reading('R1.1 kick per readout, ball layout, all spectators', lambda c: kick(c, None), 'deg')
    sgn = np.sign(np.nanmedian(kick((0, 1), None))) or 1.0
    print('  sign of the kick (carried, never folded): %+d; size = kick x sign' % sgn)
    groups = (('d<=2', d <= 2), ('d=3-4', (d >= 3) & (d <= 4)), ('d>=5', d >= 5))
    for lab, sel in groups: reading('R1.2 kick size, ball layout, %s' % lab, lambda c, sel=sel: sgn * kick(c, None)[sel], 'deg')
    ks = sgn * kick((0, 1), None)
    m, lo, hi, s = boot2(ks[d >= 5], ks[d <= 2])
    print('  %-58s %+.4g [%+.4g, %+.4g] deg -> %s' % ('R1.3 size(d>=5) - size(d<=2) [ABOVE larger far, BELOW larger near]', m, lo, hi, s))
    ok = np.isfinite(ks); bs = []
    for _ in range(B):
        j = RNG.integers(0, ok.sum(), ok.sum()); x = inv[ok][j]; y = ks[ok][j]; bs.append(np.polyfit(x, y, 1)[0])
    sl = np.polyfit(inv[ok], ks[ok], 1)[0]; lo, hi = np.percentile(bs, [2.5, 97.5])
    print('  %-58s %+.4g [%+.4g, %+.4g] deg per unit -> %s' % ('R1.4 slope of size on sum over measured of 1/d^2', sl, lo, hi, outcome(lo, hi, ok.sum())))
    d2, _, kick2 = out['spread']
    a = sgn * kick2((0, 1), None)[d2 <= 2]; b_ = ks[d <= 2]; m, lo, hi, s = boot2(a, b_)
    print('  %-58s %+.4g [%+.4g, %+.4g] deg -> %s' % ('R1.5 size(d<=2) spread layout - ball layout (same job)', m, lo, hi, s))
# ------------------------------------------------------------------------------------------------ run 2: readout count
def count(J):
    I = J.info; S = I['S']; nS = len(S); Ns = sorted({m['N'] for m in J.meta if m['N'] > 0})
    def ph(c, mode, N): return np.arctan2(J.mean_pm(c, 'fin', nS, mode=mode, N=N, setting='YY'), J.mean_pm(c, 'fin', nS, mode=mode, N=N, setting='XX'))
    def co(c, mode, N): return np.hypot(J.mean_pm(c, 'fin', nS, mode=mode, N=N, setting='YY'), J.mean_pm(c, 'fin', nS, mode=mode, N=N, setting='XX'))
    def total(c):   # phase gap M - D at each N, unwrapped along N by adding wrapped increments
        dphi = {N: wrap(ph(c, 'M', N) - ph(c, 'D', N)) for N in Ns}; tot = {Ns[0]: dphi[Ns[0]]}
        for a, b_ in zip(Ns[:-1], Ns[1:]): tot[b_] = tot[a] + wrap(dphi[b_] - dphi[a])
        return {N: np.degrees(v) for N, v in tot.items()}
    for N in Ns: reading('R2.1 kick per readout at N=%d' % N, lambda c, N=N: total(c)[N] / N, 'deg')
    reading('R2.2 kick per readout, readouts 27->81 minus readouts 1->9', lambda c: (total(c)[81] - total(c)[27]) / 54 - (total(c)[9] - total(c)[1]) / 8, 'deg')
    def icpt(c):
        t = total(c); A = np.stack([np.array(Ns, float), np.ones(len(Ns))], 1); return np.linalg.lstsq(A, np.stack([t[N] for N in Ns], 1).T, rcond=None)[0][1]
    reading('R2.3 intercept of the phase gap against N', icpt, 'deg')
    for N in Ns: reading('R2.4 coherence loss at N=%d (idle - readout)' % N, lambda c, N=N: co(c, 'D', N) - co(c, 'M', N))
    reading('R2.5 SUM form: loss(81) - 3 x loss(27)', lambda c: (co(c, 'D', 81) - co(c, 'M', 81)) - 3 * (co(c, 'D', 27) - co(c, 'M', 27)))
    def lpr(c, N):
        with np.errstate(divide='ignore', invalid='ignore'): r = co(c, 'M', N) / co(c, 'D', N); return np.where((r > 0) & (co(c, 'D', N) > 0.05), -np.log(r) / N, np.nan)
    for N in (9, 27, 81): reading('R2.6 PRODUCT form: -ln(readout/idle)/N at N=%d' % N, lambda c, N=N: lpr(c, N))
    reading('R2.7 PRODUCT form: per-readout loss at 81 minus at 9', lambda c: lpr(c, 81) - lpr(c, 9))
# ------------------------------------------------------------------------------------------------ run 3: read order
def order(J):
    G = J.info['groups']; n = sum(len(v) for v in G.values()); grp = {q: g for g, v in G.items() for q in v}
    def p_order(c, k):   # per qubit, P(1) when its group sits at slot k in the full orders
        out = np.full(n, np.nan)
        for g in 'ABC':
            perms = [o for o in ('ABC', 'ACB', 'BAC', 'BCA', 'CAB', 'CBA') if o[k] == g]
            x = np.concatenate([np.concatenate([J.get(i, 'c', n) for i in J.idx(c, kind='order', order=o)]) for o in perms])
            q = np.array(G[g]); out[q] = x[:, q].mean(0)
        return out
    def p_alone(c, k):
        out = np.full(n, np.nan)
        for g in 'ABC':
            x = np.concatenate([J.get(i, 'c', n) for i in J.idx(c, kind='alone', group=g, slot=k)]); q = np.array(G[g]); out[q] = x[:, q].mean(0)
        return out
    def p_tog(c): return np.concatenate([J.get(i, 'c', n) for i in J.idx(c, kind='together')]).mean(0)
    for k in range(3): reading('R3.1 P(1) read at slot %d after others, minus read alone at slot %d' % (k + 1, k + 1), lambda c, k=k: p_order(c, k) - p_alone(c, k))
    reading('R3.2 that excess at slot 3 minus at slot 2', lambda c: (p_order(c, 2) - p_alone(c, 2)) - (p_order(c, 1) - p_alone(c, 1)))
    reading('R3.3 read alone at slot 3 minus alone at slot 1 (waiting longer)', lambda c: p_alone(c, 2) - p_alone(c, 0))
    reading('R3.4 all read together at slot 1 minus alone at slot 1', lambda c: p_tog(c) - p_alone(c, 0))
    reading('R3.5 read last (slot 3) minus read first (slot 1), in the full orders', lambda c: p_order(c, 2) - p_order(c, 0))
# ------------------------------------------------------------------------------------------------ run 4: gates
def gates(J):
    n = J.props_n
    def P(c, ng, reps): return J.p1(c, 'c', n, arm='G', ngates=ng, reps=reps)
    for ng in (2, 6, 10): reading('R4.1 P(1) X X X train minus X train, %d steps (%d extra gates)' % (ng, 2 * ng), lambda c, ng=ng: P(c, ng, 3) - P(c, ng, 1))
    def slope(c):
        x = np.array([4., 12., 20.]); Y = np.stack([P(c, ng, 3) - P(c, ng, 1) for ng in (2, 6, 10)], 1); A = np.stack([x, np.ones(3)], 1)
        return np.linalg.lstsq(A, Y.T, rcond=None)[0][0]
    stated = np.full(n, np.nan)
    for g in J.props['gates']:
        if g['gate'] == 'x' and len(g['qubits']) == 1:
            for p in g['parameters']:
                if p['name'] == 'gate_error': stated[g['qubits'][0]] = p['value']
    reading('R4.2 cost per extra X gate (slope; 20% trimmed mean)', slope, stat='trim')
    reading('R4.3 cost per extra X gate minus IBM stated X error (20% trimmed mean)', lambda c: slope(c) - stated, stat='trim')
    print('  IBM stated X error, median over qubits: %.2e' % np.nanmedian(stated))
    for w in (30, 90, 300): reading('R4.4 reference: X then wait %d us, P(1)' % w, lambda c, w=w: J.p1(c, 'c', n, arm='G', ref_wait_us=w))
# ------------------------------------------------------------------------------------------------ run 5: neighbour state
def neighbour(J):
    T = J.info['targets']; nT = len(T); Ms = sorted({m['m'] for m in J.meta})
    def s(c, arm):
        Y = np.stack([J.p1(c, 't', nT, arm=arm, m=m) for m in Ms], 1); A = np.stack([np.array(Ms, float), np.ones(len(Ms))], 1)
        return np.linalg.lstsq(A, Y.T, rcond=None)[0][0]
    for arm in ('0', '1', '+'): reading('R5.1 error per X gate, neighbours in %s' % arm, lambda c, arm=arm: s(c, arm))
    reading('R5.2 neighbours in 1 minus in 0 (static neighbour state)', lambda c: s(c, '1') - s(c, '0'))
    reading('R5.3 neighbours in + minus in 0 (the original comparison)', lambda c: s(c, '+') - s(c, '0'))
    reading('R5.4 neighbours in + minus the mean of 0 and 1', lambda c: s(c, '+') - (s(c, '0') + s(c, '1')) / 2)
    with np.errstate(divide='ignore', invalid='ignore'): r = s((0, 1), '+') / s((0, 1), '0')
    m, lo, hi, k, o = boot(r); print('  R5.5 ratio + / 0 (beside, the original 1.7x): %.3f [%.3f, %.3f] n=%d' % (m, lo, hi, k))
# ------------------------------------------------------------------------------------------------ run 6: revival
DIP, REV = (26, 50), (58, 78)
def fit_heldout(t, y):
    """Per qubit: fit on odd-index taus, score on even-index taus. Plain decay A e^(-t/T) + c vs damped beat A e^(-t/T)|cos(pi t/P)| + c.
    Grid over T (and P); A, c by least squares. Returns held-out SSE(decay) - SSE(beat)."""
    tr = np.arange(len(t)) % 2 == 1; te = ~tr; Ts = np.geomspace(5, 500, 40); Ps = np.arange(40, 101, 1.0)
    def best(shapes):
        bs, bp = np.inf, None
        for f in shapes:
            A = np.stack([f[tr], np.ones(tr.sum())], 1); co, *_ = np.linalg.lstsq(A, y[tr], rcond=None); r = np.sum((A @ co - y[tr]) ** 2)
            if r < bs: bs, bp = r, (f, co)
        f, co = bp; return np.sum((co[0] * f[te] + co[1] - y[te]) ** 2)
    dec = best([np.exp(-t / T) for T in Ts]); beat = best([np.exp(-t / T) * np.abs(np.cos(np.pi * t / P)) for T in Ts for P in Ps])
    return dec - beat
def revival(J):
    A = J.info['A']; nA = len(A); taus = np.array(sorted({m['tau'] for m in J.meta}), float)
    def C(c, arm): return np.stack([np.hypot(J.mean_pm(c, 'a', nA, arm=arm, tau=int(t), setting='X'), J.mean_pm(c, 'a', nA, arm=arm, tau=int(t), setting='Y')) for t in taus], 1)
    def R(c, arm):
        x = C(c, arm); w = lambda lo, hi: x[:, (taus >= lo) & (taus <= hi)].mean(1); return w(*REV) - w(*DIP)
    for arm in ('0', '1', '+'):
        reading('R6.1 come-back: mean C(58-78 us) - mean C(26-50 us), B in %s' % arm, lambda c, arm=arm: R(c, arm))
    reading('R6.2 come-back, B in + minus B in 0', lambda c: R(c, '+') - R(c, '0'))
    reading('R6.3 come-back, B in 1 minus B in 0', lambda c: R(c, '1') - R(c, '0'))
    for arm in ('0', '1', '+'):
        reading('R6.4 held-out SSE(plain decay) - SSE(damped beat), B in %s' % arm, lambda c, arm=arm: np.array([fit_heldout(taus, y) for y in C(c, arm)]))
def summary():
    """Chips combined by VALUE and WIDTH, never by label: per chip value [interval]; inverse-variance combination (sigma from the 95 % interval
    width / 3.92) with its interval and outcome; and the chip-to-chip spread, chi2 of the chips about the combination over its dof, as a reading
    of its own (about 1 when the chips agree within their widths, larger when they differ)."""
    print('=' * 100); print('ACROSS CHIPS (values and widths; a chip whose copies shifted is listed and left out of the combination)')
    by = {}
    for lab, chip, p in LOG: by.setdefault(lab, []).append((chip, p))
    for lab, rows in by.items():
        ok = [(c, p) for c, p in rows if p is not None and np.isfinite(p[1])]
        cells = ' | '.join('%s %+.4g [%+.4g, %+.4g]' % (c.replace('ibm_', ''), p[0], p[1], p[2]) if p is not None else '%s copies shifted' % c.replace('ibm_', '') for c, p in rows)
        if len(ok) < 2: print('  %s\n      %s' % (lab, cells)); continue
        v = np.array([p[0] for _, p in ok]); s = np.array([(p[2] - p[1]) / 3.92 for _, p in ok]); s = np.where(s > 0, s, np.nan)
        w = 1 / s ** 2; m = np.nansum(w * v) / np.nansum(w); se = 1 / np.sqrt(np.nansum(w)); lo, hi = m - 1.96 * se, m + 1.96 * se
        chi = np.nansum(((v - m) / s) ** 2) / (len(ok) - 1)
        print('  %s\n      %s\n      combined %+.4g [%+.4g, %+.4g] -> %s | chip-to-chip spread chi2/dof %.2f (dof %d)' % (lab, cells, m, lo, hi, outcome(lo, hi, 99), chi, len(ok) - 1))
RUNS = dict(distance=distance, count=count, order=order, gates=gates, neighbour=neighbour, revival=revival)
def run_all(folder=HERE, synth=None):
    for jf in sorted(glob.glob(os.path.join(folder, 'job_*_*.json'))):
        run, chip = os.path.basename(jf)[4:-5].split('_', 1)
        if not os.path.exists(os.path.join(folder, 'raw_%s_%s.npz' % (run, chip))) and synth is None: print('=' * 100); print(run, chip, 'NOTHING RETURNED'); continue
        J = Job(run, chip, folder, None if synth is None else synth[(run, chip)]); J.props_n = J.rec.get('n_qubits', 156)
        print('=' * 100); print(run, chip, 'job', J.rec['job_id']); CUR[0] = chip; RUNS[run](J)
    summary()
if __name__ == '__main__': run_all()
