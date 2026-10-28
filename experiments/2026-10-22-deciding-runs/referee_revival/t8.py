from model2 import *
import model, model2
rng=np.random.default_rng(2)
rec,props,A,B,taus,M=load('phoenix'); E_,zz=zzmap(props)
zz2={k:v*0.5 for k,v in zz.items()}; model2.zzmap=lambda p:(E_,zz2); model.zzmap=model2.zzmap
res,_,nb,_=model2.run('phoenix',verbose=False)
Rp,Rm,pred,C=res['+']; deg=np.array([len(nb[a]) for a in A])
def tc(Cc): return np.array([taus[np.argmax(c<0.2)] for c in Cc])
for d in (3,4):
    m=deg==d
    print('deg',d,'n',m.sum(),'meas: t(C+<0.2) median %.0f us, R+ %+.4f | pred(scale .5): t %.0f R %+.4f'%(np.median(tc(C[m])),np.median(Rm[m]),np.median(tc(pred[m])),np.median(Rp[m])))
m3=deg==3;m4=deg==4
b=[np.median(rng.choice(tc(C[m4]),m4.sum()))-np.median(rng.choice(tc(C[m3]),m3.sum())) for _ in range(4000)]
print('t(C<0.2) deg4-deg3 meas %.1f [%.1f,%.1f]'%(np.median(tc(C[m4]))-np.median(tc(C[m3])),*np.percentile(b,[2.5,97.5])))
# bump specific: peak of C+ in 54-82 minus min in 30-54, per qubit
w=lambda Z,lo,hi:Z[:,(taus>=lo)&(taus<=hi)]
bump=w(C,56,82).max(1)-w(C,30,54).min(1); print('per-qubit bump (max C 56-82 minus min C 30-54): median %.3f; same stat on arm 0 and 1 not meaningful'%np.median(bump))
# null of that stat from synthetic pure noise floor: C=hypot noise only
sim=np.hypot(rng.normal(0,0.0316,(60,38)),rng.normal(0,0.0316,(60,38))); 
print(' pure-noise-floor same stat median %.3f'%np.median(w(sim,56,82).max(1)-w(sim,30,54).min(1)))
# R+ from tail-only floor: R of pure floor = 0 by construction. Fraction of qubits with R+>0:
print('fraction of qubits with R+>0: %.2f (n=%d); sign test p~'%((Rm>0).mean(),len(Rm)))
from scipy.stats import binomtest; print(binomtest(int((Rm>0).sum()),len(Rm)).pvalue)
