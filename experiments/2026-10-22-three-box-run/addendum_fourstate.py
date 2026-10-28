"""addendum_fourstate.py -- ADDENDUM to PLAN_THREEBOX.md, written after the fez job was submitted and BEFORE any result was read.
Why: the plan asked 'sign OR corner' -- two lanes. The shortfall can follow the sign, the corner, BOTH, or NEITHER (another feature, e.g.
the |00> position, which needs the most X gates). This reading carries all of them.
Per chip, over the 9 support boxes (3 variants x 3 boxes): r_b = hardware P / model P, regressed as
    r_b = c0 + a*[box carries the negative sign] + b*[box is |11>] + c*[box is |00>]
a = SIGN term, b = CORNER term, c = |00>-position term, 95% intervals from the per-box bootstrap draws (2000).
State per chip from a and b only (c is reported, not used to name the state):
    SIGN (a below 0, b consistent) | CORNER (b below, a consistent) | BOTH (a and b below) | NEITHER (a and b consistent) |
    CANNOT_TELL (the sign or corner term differs between copy 0 and copy 1: copy-gap interval excludes 0) | any ABOVE reported as it falls."""
import sys, json, numpy as np
import analyze_threebox as A
from tb_common import VARIANTS, EPS, xp_from_bits, fit_SP
def _coef(X, pts): return np.linalg.lstsq(X, pts, rcond=None)[0]
def fourstate(R, raw=None):
    """R from analyze_threebox.analyze; raw = the same npz path (for the per-copy terms)."""
    rows, draws, cells = [], [], []
    for v, V in VARIANTS.items():
        for b in V['support']:
            r = R['box'][(v, b)]; rows.append([1.0, float(V['sign'][b] < 0), float(b == '11'), float(b == '00')])
            draws.append(r['D'] / r['Pm']); cells.append((v, b, r['Pm']))
    X = np.array(rows); Y = np.array(draws).T
    pt = np.array([R['box'][(v, b)]['P'] / R['box'][(v, b)]['Pm'] for v, V in VARIANTS.items() for b in V['support']])
    coef = _coef(X, pt); cd = np.array([_coef(X, y) for y in Y]); ci = np.percentile(cd, [2.5, 97.5], axis=0)
    gaps = None
    if raw is not None:
        z = np.load(raw); per = []
        for cp in (0, 1):
            dd = []
            for v, b, Pm in cells:
                arrs = [xp_from_bits(z[f'{v}|{b}|{e:.2f}|{cp}']) for e in EPS]; _, D, *_ = A.boot_P(arrs); dd.append(D / Pm)
            per.append(np.array([_coef(X, y) for y in np.array(dd).T]))
        gaps = np.percentile(per[0] - per[1], [2.5, 97.5], axis=0)
    lab = ['ABOVE' if ci[0, j] > 0 else 'BELOW' if ci[1, j] < 0 else 'CONSISTENT' for j in range(4)]
    bad = gaps is not None and any(gaps[0, j] > 0 or gaps[1, j] < 0 for j in (1, 2))
    a, b = lab[1], lab[2]
    state = ('CANNOT_TELL (copies disagree on the sign or corner term)' if bad else 'SIGN' if (a, b) == ('BELOW', 'CONSISTENT') else 'CORNER' if (a, b) == ('CONSISTENT', 'BELOW')
             else 'BOTH' if (a, b) == ('BELOW', 'BELOW') else 'NEITHER' if (a, b) == ('CONSISTENT', 'CONSISTENT') else f'OTHER: sign {a}, corner {b}')
    names = ['c0', 'a SIGN', 'b CORNER', 'c |00>']
    g = '' if gaps is None else '  copy gaps: sign [%+.3f, %+.3f], corner [%+.3f, %+.3f]' % (gaps[0, 1], gaps[1, 1], gaps[0, 2], gaps[1, 2])
    return '\n'.join([f"  {R['chip']}: " + ' | '.join(f'{n} {coef[j]:+.3f} [{ci[0,j]:+.3f}, {ci[1,j]:+.3f}] {lab[j]}' for j, n in enumerate(names)), f'  {R["chip"]} STATE: {state}{g}'])
if __name__ == '__main__':
    a = sys.argv[1:]; res = [A.analyze(a[i], json.load(open(a[i + 1]))) for i in range(0, len(a), 2)]
    print('FOUR-STATE READING (addendum)'); [print(fourstate(R, a[2 * i])) for i, R in enumerate(res)]
    print('R4 common probe offset: mean <Xp> on kept shots of the unused box over all eps, hardware vs model')
    for i, R in enumerate(res):
        z = np.load(a[2 * i]); m = json.load(open(a[2 * i + 1]))['model']
        for v, V in VARIANTS.items():
            ub = [b for b in ('00', '01', '10', '11') if b not in V['support']][0]
            x = np.concatenate([xp_from_bits(z[f'{v}|{ub}|{e:.2f}|{cp}']) for e in EPS for cp in (0, 1)])
            mm = float(np.mean([m[f'{v}|{ub}|{e:.2f}'] for e in EPS])); se = x.std(ddof=1) / np.sqrt(len(x))
            print(f"  {R['chip']} {v} unused {ub}: hardware {x.mean():+.4f} +- {se:.4f} (n {len(x)}), model {mm:+.4f}, difference {x.mean() - mm:+.4f} [{x.mean() - mm - 1.96*se:+.4f}, {x.mean() - mm + 1.96*se:+.4f}]")
