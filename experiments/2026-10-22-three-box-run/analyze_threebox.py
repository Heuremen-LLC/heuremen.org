"""analyze_threebox.py -- committed analysis for the three-box run (PLAN_THREEBOX.md). Usage: analyze_threebox.py <raw.npz> <model.json> [<raw> <model> ...]
raw npz: arrays '<var>|<box>|<eps>|<copy>' of packed bits (bit0 s0, bit1 s1, bit2 probe). model json: {'chip':..,'model':{key:<Xp>},'ideal':{key:<Xp>}}"""
import sys, json, math, numpy as np
from tb_common import EPS, BOXES, VARIANTS, textbook, xp_from_bits, fit_SP
NB = 2000; MINKEPT = 50; rng = np.random.default_rng(20261022)
def k(v, b, e): return f'{v}|{b}|{e:.2f}'
def outcome(lo, hi, gap=None, enough=True):
    if not enough: return 'CANNOT_TELL (too few kept shots)'
    if gap is not None and (gap[0] > 0 or gap[1] < 0): return 'CANNOT_TELL (copies disagree)'
    return 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def boot_P(arrs):
    """arrs: list over EPS of +-1 arrays. Returns P point, P bootstrap draws (binomial resampling within each circuit), S point, chi2."""
    n = np.array([len(a) for a in arrs]); m = np.array([a.mean() if len(a) else 0 for a in arrs]); se = np.sqrt(np.maximum(1 - m ** 2, 1e-6) / np.maximum(n, 1))
    S, P, chi = fit_SP(EPS, m, se); ph = (1 + m) / 2; draws = []
    for _ in range(NB):
        mb = 2 * rng.binomial(n, ph) / np.maximum(n, 1) - 1; draws.append(fit_SP(EPS, mb, se)[1])
    return P, np.array(draws), S, chi, int(n.min())
def analyze(raw, model):
    z = np.load(raw); chip = model['chip']; out = {'chip': chip, 'box': {}}
    for v in VARIANTS:
        for b in BOXES:
            c = [[xp_from_bits(z[f'{k(v,b,e)}|{cp}']) for e in EPS] for cp in (0, 1)]
            P, D, S, chi, nmin = boot_P([np.concatenate([c[0][i], c[1][i]]) for i in range(len(EPS))])
            P0, D0, *_ = boot_P(c[0]); P1, D1, *_ = boot_P(c[1])
            gap = np.percentile(D0 - D1, [2.5, 97.5])
            mm = [model['model'][k(v, b, e)] for e in EPS]; ii = [model['ideal'][k(v, b, e)] for e in EPS]
            Pm = fit_SP(EPS, mm, np.ones(len(EPS)))[1]; Pi = fit_SP(EPS, ii, np.ones(len(EPS)))[1]
            out['box'][(v, b)] = dict(P=P, D=D, S=S, chi=chi, nmin=nmin, gap=gap, Pm=Pm, Pi=Pi, tb=textbook(v)[b])
    return out
def report(results):
    L = []
    for R in results:
        L.append(f"\n=== {R['chip']} ===\nvariant box  textbook  ideal(fit)  IBM-model(fit)  HARDWARE P [95%]            S      chi2/dof  copy gap          kept/circuit")
        for (v, b), r in R['box'].items():
            lo, hi = np.percentile(r['D'], [2.5, 97.5])
            L.append(f"{v:4s} {b}   {r['tb']:+.0f}       {r['Pi']:+.3f}      {r['Pm']:+.3f}         {r['P']:+.3f} [{lo:+.3f}, {hi:+.3f}]   {r['S']:+.3f}  {r['chi']:.2f}     [{r['gap'][0]:+.3f}, {r['gap'][1]:+.3f}]  {r['nmin']}")
        L.append('R1  negative box, hardware P minus IBM-model P (decision reading per variant):')
        for v, V in VARIANTS.items():
            nb = [b for b, s in V['sign'].items() if s < 0][0]; r = R['box'][(v, nb)]; d = r['D'] - r['Pm']; lo, hi = np.percentile(d, [2.5, 97.5])
            L.append(f"    {v} box {nb}: {r['P'] - r['Pm']:+.3f} [{lo:+.3f}, {hi:+.3f}] -> {outcome(lo, hi, r['gap'], r['nmin'] >= MINKEPT)}")
        L.append('R2  sign or position: G = (hardware/IBM-model on the negative box) minus (mean of the same ratio on the two positive boxes):')
        for v, V in VARIANTS.items():
            nb = [b for b, s in V['sign'].items() if s < 0][0]; pb = [b for b, s in V['sign'].items() if s > 0]
            g = R['box'][(v, nb)]['D'] / R['box'][(v, nb)]['Pm'] - np.mean([R['box'][(v, b)]['D'] / R['box'][(v, b)]['Pm'] for b in pb], axis=0)
            gp = R['box'][(v, nb)]['P'] / R['box'][(v, nb)]['Pm'] - np.mean([R['box'][(v, b)]['P'] / R['box'][(v, b)]['Pm'] for b in pb])
            lo, hi = np.percentile(g, [2.5, 97.5]); gaps = [R['box'][(v, b)]['gap'] for b in [nb] + pb]
            bad = any(x[0] > 0 or x[1] < 0 for x in gaps)
            L.append(f"    {v} (negative at {nb}): {gp:+.3f} [{lo:+.3f}, {hi:+.3f}] -> {outcome(lo, hi, None, not bad) if not bad else 'CANNOT_TELL (copies disagree)'}")
        for v in VARIANTS:
            ub = [b for b in BOXES if b not in VARIANTS[v]['support']][0]; r = R['box'][(v, ub)]
            L.append(f"R3  {v} unused box {ub}: S offset {r['S']:+.3f}, P {r['P']:+.3f} (IBM-model P {r['Pm']:+.3f})")
    return '\n'.join(L)
def combined(results):
    """inverse-variance combination of R2 across chips per variant; quoted only if chi2/dof < 2."""
    L = ['\nCOMBINED across chips (R2, per variant; sigma = 95% width / 3.92):']
    for v, V in VARIANTS.items():
        nb = [b for b, s in V['sign'].items() if s < 0][0]; pb = [b for b, s in V['sign'].items() if s > 0]; vals = []
        for R in results:
            g = R['box'][(v, nb)]['D'] / R['box'][(v, nb)]['Pm'] - np.mean([R['box'][(v, b)]['D'] / R['box'][(v, b)]['Pm'] for b in pb], axis=0)
            gp = R['box'][(v, nb)]['P'] / R['box'][(v, nb)]['Pm'] - np.mean([R['box'][(v, b)]['P'] / R['box'][(v, b)]['Pm'] for b in pb])
            lo, hi = np.percentile(g, [2.5, 97.5]); vals.append((gp, (hi - lo) / 3.92))
        w = np.array([1 / s ** 2 for _, s in vals]); m = float((w * np.array([x for x, _ in vals])).sum() / w.sum()); s = float(1 / math.sqrt(w.sum()))
        chi = float((w * (np.array([x for x, _ in vals]) - m) ** 2).sum() / max(1, len(vals) - 1))
        lab = outcome(m - 1.96 * s, m + 1.96 * s) if chi < 2 else 'not combined (chip-to-chip chi2/dof >= 2)'
        L.append(f"    {v}: {m:+.3f} [{m - 1.96*s:+.3f}, {m + 1.96*s:+.3f}] chi2/dof {chi:.2f} -> {lab}")
    return '\n'.join(L)
if __name__ == '__main__':
    a = sys.argv[1:]; res = [analyze(a[i], json.load(open(a[i + 1]))) for i in range(0, len(a), 2)]
    print(report(res)); print(combined(res))
