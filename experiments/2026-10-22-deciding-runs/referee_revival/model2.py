from model import *
from scipy.optimize import curve_fit
def floor(c): return np.sqrt(c**2+0.0021)  # ~hypot noise bias for 1000 shots
def run(chip,zscale=1.0,verbose=True):
    rec,props,A,B,taus,M=load(chip)
    A,B,nb,zz,P=predict(chip,taus)
    X0,Y0,C0=curves(M,taus,'0')
    # T2* from arm0 : C=a*exp(-(t/T)^p)  p=1 via log fit on first part
    t=taus*1e-6; res={}
    # extra dephasing factor = measured arm0 / (nothing): use smoothed arm0 curve per qubit (fit exp)
    ext=[]
    for i,a in enumerate(A):
        f=lambda tt,amp,T,p: amp*np.exp(-(tt/T)**p)
        try: po,_=curve_fit(f,t,C0[i],p0=[0.97,100e-6,1.0],bounds=([0.5,5e-6,0.5],[1.2,1e-3,2.5]))
        except Exception: po=[0.97,100e-6,1.0]
        ext.append(f(t,*po))
    ext=np.array(ext)
    # divide out props T2 part we put into P: recompute P with no T2 -> P/ exp(-t/T2)
    T2=qp(props,'T2')*1e-6
    for arm in '01+':
        Pz=P[arm]/np.exp(-t[None,:]/T2[A][:,None])
        pred=floor(ext*Pz)
        X,Y,C=curves(M,taus,arm)
        w=lambda Z,lo,hi:Z[:,(taus>=lo)&(taus<=hi)].mean(1)
        Rm=w(C,58,78)-w(C,26,50); Rp=w(pred,58,78)-w(pred,26,50)
        res[arm]=(Rp,Rm,pred,C)
        if verbose: print(chip,'arm',arm,'pred med R %+.4f  meas %+.4f  per-qubit corr %.3f  rmse %.3f'%(np.median(Rp),np.median(Rm),np.corrcoef(Rp,Rm)[0,1],np.sqrt(np.mean((Rp-Rm)**2))))
        if verbose and arm=='+':
            print('   pred C+',' '.join('%5.2f'%v for v in np.median(pred,0)[10:26])); print('   meas C+',' '.join('%5.2f'%v for v in np.median(C,0)[10:26])); print('   taus  ',' '.join('%5d'%v for v in taus[10:26]))
    return res,A,nb,taus
if __name__=='__main__':
    for chip in ('phoenix','fez'): run(chip)
