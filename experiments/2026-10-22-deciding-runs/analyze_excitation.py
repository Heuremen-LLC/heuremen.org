"""Analysis exactly as PLAN_EXCITATION.md. Uses analyze.py's readers and rulers unchanged."""
import os, sys, glob
import numpy as np
import analyze as A

def excitation(J):
    nS = len(J.info['S'])
    def s(c, mo):   # N = 0 has only the idle circuit for each state, shared by D and M
        out = []
        for N in (0, 4, 8):
            m = ('D' + mo[1]) if N == 0 else mo
            out.append(np.arctan2(J.mean_pm(c, 'fin', nS, mode=m, N=N, setting='YY'), J.mean_pm(c, 'fin', nS, mode=m, N=N, setting='XX')))
        return np.degrees(np.polyfit([0., 4., 8.], np.unwrap(np.array(out), axis=0), 1)[0])
    A.reading('RJ.1 kick, measured qubits in 0 (M0 - D0)', lambda c: s(c, 'M0') - s(c, 'D0'), 'deg/step')
    A.reading('RJ.2 kick, measured qubits in 1 (M1 - D1)', lambda c: s(c, 'M1') - s(c, 'D1'), 'deg/step')
    A.reading('RJ.3 G = kick in 1 minus kick in 0', lambda c: (s(c, 'M1') - s(c, 'D1')) - (s(c, 'M0') - s(c, 'D0')), 'deg/step')
    A.reading('RJ.4 measured qubits in 1, no readout (D1 - D0)', lambda c: s(c, 'D1') - s(c, 'D0'), 'deg/step')
A.RUNS['excitation'] = excitation
if __name__ == '__main__':
    folder = sys.argv[1] if len(sys.argv) > 1 else A.HERE
    for jf in sorted(glob.glob(os.path.join(folder, 'job_excitation_*.json'))):
        chip = os.path.basename(jf)[len('job_excitation_'):-5]
        if not os.path.exists(os.path.join(folder, 'raw_excitation_%s.npz' % chip)): print(chip, 'NOTHING RETURNED'); continue
        J = A.Job('excitation', chip, folder); J.props_n = J.rec.get('n_qubits', 156); print('=' * 100); print('excitation', chip, 'job', J.rec['job_id']); A.CUR[0] = chip; excitation(J)
    A.summary()
