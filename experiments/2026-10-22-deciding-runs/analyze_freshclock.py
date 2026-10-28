"""Analysis exactly as PLAN_FRESHCLOCK.md. Uses analyze.py's readers and rulers unchanged."""
import os, sys, glob
import numpy as np
import analyze as A
def freshclock(J):
    nS = len(J.info['S']); W = J.info['window']
    ph = lambda c, pre, win: np.arctan2(J.mean_pm(c, 'fin', nS, pre=pre, win=win, setting='YY'), J.mean_pm(c, 'fin', nS, pre=pre, win=win, setting='XX'))
    co = lambda c, pre, win: np.hypot(J.mean_pm(c, 'fin', nS, pre=pre, win=win, setting='YY'), J.mean_pm(c, 'fin', nS, pre=pre, win=win, setting='XX'))
    kick = lambda c, pre: np.degrees(A.wrap(ph(c, pre, 'M') - ph(c, pre, 'D'))) / W
    for pre in ('M0', 'M27', 'M72', 'D72'):
        print('  instrument: window idle coherence after preload %s, median %.3f' % (pre, np.median(co((0, 1), pre, 'D'))))
    sgn = np.sign(np.nanmedian(kick((0, 1), 'M0'))) or 1.0
    print('  sign of the kick after no preload: %+d; size = kick x sign (carried, never folded)' % sgn)
    for pre, lab in (('M0', 'no preload'), ('M27', 'after 27 readouts'), ('M72', 'after 72 readouts'), ('D72', 'after 72 idles')):
        A.reading('RF.1 kick size per readout, fresh clock, %s' % lab, lambda c, pre=pre: sgn * kick(c, pre), 'deg')
    A.reading('RF.2 size after 72 readouts minus after none (BELOW: fades)', lambda c: sgn * (kick(c, 'M72') - kick(c, 'M0')), 'deg')
    A.reading('RF.3 size after 27 readouts minus after none', lambda c: sgn * (kick(c, 'M27') - kick(c, 'M0')), 'deg')
    A.reading('RF.4 size after 72 readouts minus after 72 idles (readouts, not time)', lambda c: sgn * (kick(c, 'M72') - kick(c, 'D72')), 'deg')
A.RUNS['freshclock'] = freshclock
if __name__ == '__main__':
    folder = sys.argv[1] if len(sys.argv) > 1 else A.HERE
    for jf in sorted(glob.glob(os.path.join(folder, 'job_freshclock_*.json'))):
        chip = os.path.basename(jf)[len('job_freshclock_'):-5]
        if not os.path.exists(os.path.join(folder, 'raw_freshclock_%s.npz' % chip)): print(chip, 'NOTHING RETURNED'); continue
        J = A.Job('freshclock', chip, folder); J.props_n = J.rec.get('n_qubits', 156); print('=' * 100); print('freshclock', chip, 'job', J.rec['job_id']); A.CUR[0] = chip; freshclock(J)
    A.summary()
