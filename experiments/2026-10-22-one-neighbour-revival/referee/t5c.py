import sys;sys.path.insert(0,'/private/tmp/claude-501/-Users-heuremen/80087010-0cb3-4ae5-9e13-362e6d5086ad/scratchpad/ref_zz')
from lib import *
G=np.arange(0.10,2.5001,0.0025)
rng=np.random.default_rng(5)
for chip in('ibm_phoenix','ibm_fez'):
    Jb=J(chip);ts=np.array(Jb.taus)*1e-6
    C={a:np.array([coh(Jb,a,t,(0,1)) for t in Jb.taus]) for a in('Z','P1','P2')}
    R=[]
    for j,(b1,a,b2) in enumerate(Jb.tr):
        for arm,q in(('P1',b1),('P2',b2)):
            f=Jb.zz[tuple(sorted((a,q)))]
            if f<=0:continue
            mod=C['Z'][:,j][None]*np.abs(np.cos(np.pi*G[:,None]*f*ts[None]));k=np.argmin(((mod-C[arm][:,j][None])**2).sum(1))
            # free-frequency fit too: f grid 1-100 kHz
            F=np.arange(1e3,100e3,250.);mod=C['Z'][:,j][None]*np.abs(np.cos(np.pi*F[:,None]*ts[None]));kf=np.argmin(((mod-C[arm][:,j][None])**2).sum(1))
            sse=((mod[kf]-C[arm][:,j])**2).sum()
            R.append((j,f/1e3,G[k]*f/1e3,F[kf]/1e3,sse))
    R=np.array(R);ok=R[:,3]<99;lf=np.log(R[ok,1]);lb=np.log(R[ok,3])
    sl=np.polyfit(lf,lb,1)[0];bs=[]
    jj=np.where(ok)[0]
    for _ in range(2000):
        i=rng.choice(len(jj),len(jj));bs.append(np.polyfit(lf[i],lb[i],1)[0])
    print(chip,'n pairs',ok.sum(),'log-log slope of free fitted beat vs props zz = %.2f [%.2f, %.2f] (1=tracks, 0=constant)'%(sl,*np.percentile(bs,[2.5,97.5])))
    print('  free-fit beat kHz: median %.1f IQR %.1f-%.1f ; props zz median %.1f IQR %.1f-%.1f ; corr %.2f'%(np.median(R[ok,3]),*np.percentile(R[ok,3],[25,75]),np.median(R[ok,1]),*np.percentile(R[ok,1],[25,75]),np.corrcoef(R[ok,1],R[ok,3])[0,1]))
    # per-pair ratio free/zz median
    print('  median free/zz = %.2f'%np.median(R[ok,3]/R[ok,1]))
    rk=lambda x:np.argsort(np.argsort(x))
    print('  Spearman free-fit beat vs props zz: %.2f (n=%d)'%(np.corrcoef(rk(R[ok,1]),rk(R[ok,3]))[0,1],ok.sum()))
