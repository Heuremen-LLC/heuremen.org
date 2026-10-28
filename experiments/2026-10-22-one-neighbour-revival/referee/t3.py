import sys;sys.path.insert(0,'/Users/heuremen/bones-forge/arxiv_one_readout/zz_revival')
import analyze_zz as A, numpy as np
rng=np.random.default_rng(3)
def med(x):
    x=np.asarray(x,float);x=x[np.isfinite(x)]
    if len(x)<5:return 'n=%d'%len(x)
    b=np.median(rng.choice(x,(2000,len(x))),1);return '%+.3f [%+.3f, %+.3f] n=%d'%(np.median(x),*np.percentile(b,[2.5,97.5]),len(x))
for chip in('ibm_fez','ibm_phoenix'):
    Jb=A.Job(chip);T=np.array(Jb.taus);m=Jb.m
    C={a:np.array([A.coh(Jb,a,t,(0,1))[0] for t in Jb.taus]) for a in('Z','P1')}
    zz=A.zz_table(Jb.props,Jb.edges);f1=np.array([zz[tuple(sorted((a,b1)))] for b1,a,b2 in Jb.trip])
    R,info=A.readings(Jb,(0,1))
    h=R['R2 return height at revival / Z (as measured)'];d=R['R2 return height at revival / Z, minus ZZ+T1 prediction (1+exp(-t*/T1_B1))/2']
    s=R['R1 fitted scale s (P1, B1 alone in +)']
    T1=np.array([A.qp(Jb.props,b1,'T1') for b1,a,b2 in Jb.trip])
    ts=1/(s*f1)*1e6
    floor=5*np.sqrt(2/(np.pi*2*1000))
    zc=np.full(m,np.nan);pc=np.full(m,np.nan)
    for j in range(m):
        if np.isfinite(h[j]):
            k=np.argsort(abs(T-ts[j]))[:3];zc[j]=C['Z'][k,j].mean();pc[j]=C['P1'][k,j].mean()
    print('=====',chip,'noise floor(5x)=%.3f  1x=%.3f  (rice floor of C ~ %.3f)'%(floor,floor/5,np.sqrt(np.pi/ (4*2000))*2))
    print(' R2 as committed',med(d))
    ok=np.isfinite(d)
    print(' t* within last 10%% of grid (>%.0f us): n=%d of %d'%(T[0]+0.9*(T[-1]-T[0]),(ok&(ts>T[0]+0.9*(T[-1]-T[0]))).sum(),ok.sum()))
    for name,mask in (('drop t*>273us',ts<=T[0]+0.9*(T[-1]-T[0])),('Z(t*)>0.15',zc>0.15),('Z(t*)>0.25',zc>0.25),('Z(t*)>0.35',zc>0.35),('Z(t*)<=0.15',zc<=0.15),('t*<=T[-2]-keep & fitted s in 0.5-1.5',(s>0.5)&(s<1.5)),('t*<150us',ts<150)):
        print('  %-40s R2 diff %s | h %s'%(name,med(np.where(mask&ok,d,np.nan)),med(np.where(mask&ok,h,np.nan))))
    print(' per-triple (fez): t*  Z(t*)  P1(t*)  h  hp  d  T1_B1')
    for j in np.argsort(ts):
        if ok[j]:print('   %2d t*=%6.1f Z=%.2f P1=%.2f h=%.2f pred=%.2f d=%+.2f T1B1=%5.1f s=%.2f'%(j,ts[j],zc[j],pc[j],h[j],h[j]-d[j],d[j],T1[j],s[j]))
