import sys;sys.path.insert(0,'/Users/heuremen/bones-forge/arxiv_one_readout/zz_revival')
import analyze_zz as A, numpy as np
from scipy.optimize import curve_fit
for chip in('ibm_fez','ibm_phoenix'):
    Jb=A.Job(chip);T=np.array(Jb.taus);m=Jb.m
    # P(b1=1 | O1) vs tau, both copies, X and Y settings pooled
    P=np.zeros((len(T),m))
    for ti,t in enumerate(Jb.taus):
        acc=[];
        for c in(0,1):
            for st in 'XY':
                i=Jb.idx(arm='O1',tau_us=t,setting=st,copy=c)[0];acc.append(Jb.reg('b1',i).mean(0))
        P[ti]=np.mean(acc,0)
    Z0=np.zeros((len(T),m))
    for ti,t in enumerate(Jb.taus):
        Z0[ti]=np.mean([Jb.reg('b1',Jb.idx(arm='Z',tau_us=t,setting=st,copy=c)[0]).mean(0) for c in(0,1) for st in 'XY'],0)
    T1p=np.array([A.qp(Jb.props,b1,'T1') for b1,a,b2 in Jb.trip]);T1m=np.full(m,np.nan)
    for j in range(m):
        try:
            p,_=curve_fit(lambda t,a,T1,c:a*np.exp(-t/T1)+c,T,P[:,j],p0=[0.95,100,0.02],bounds=([0.3,5,0],[1.01,2000,0.2]));T1m[j]=p[1]
        except Exception as e:pass
    print(chip,'B1 T1 props median %.0f us, measured from O1 b1 decay median %.0f us; ratio meas/props median %.2f (IQR %.2f-%.2f); P(b1=1) in Z arm, tau-mean per triple median %.3f max %.3f'%(np.median(T1p),np.nanmedian(T1m),np.nanmedian(T1m/T1p),*np.nanpercentile(T1m/T1p,[25,75]),np.median(Z0.mean(0)),Z0.mean(0).max()))
    R,info=A.readings(Jb,(0,1))
    h=R['R2 return height at revival / Z (as measured)'];s=R['R1 fitted scale s (P1, B1 alone in +)'];ts=1/(s*info['f1'])*1e6
    hp_m=(1+np.exp(-ts/T1m))/2
    d=h-hp_m;print('  R2 with measured T1 of B1: median %.3f ; n=%d;  with T1=inf bound (h-1): %.3f'%(np.nanmedian(d),np.isfinite(d).sum(),np.nanmedian(h-1)))
    # Z-arm excited population of B2 as well
