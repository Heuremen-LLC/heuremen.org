from model2 import *
rng=np.random.default_rng(3)
rec,props,A,B,taus,M=load('phoenix'); meta=rec['meta']
# 3: copy shift in arm 1: by tau window and by execution order half
def Rof(arm,pred):
    # pred(meta_i)->bool, rebuild means from raw
    z=np.load(D+'raw_revival_ibm_phoenix.npz'); acc={}
    for i,m in enumerate(meta):
        if m['arm']!=arm or not pred(i,m): continue
        acc.setdefault((m['tau'],m['setting']),[]).append(1-2*bits(z['a_%d'%i],len(A)))
    Cc=np.stack([np.hypot(np.concatenate(acc[(int(t),'X')]).mean(0),np.concatenate(acc[(int(t),'Y')]).mean(0)) if (int(t),'X') in acc and (int(t),'Y') in acc else np.full(len(A),np.nan) for t in taus],1)
    w=lambda lo,hi:np.nanmean(Cc[:,(taus>=lo)&(taus<=hi)],1); return w(58,78)-w(26,50),Cc
order=rec['order']; pos=np.empty(len(order),int)
# meta index i is position in the shuffled pub list = execution order
n=len(meta)
for arm in '01+':
    for lab,pred in (('copy0',lambda i,m:m['copy']==0),('copy1',lambda i,m:m['copy']==1),('first half exec',lambda i,m:i<n//2),('second half exec',lambda i,m:i>=n//2)):
        r,_=Rof(arm,pred); print('arm',arm,lab,'median R %+.4f'%np.median(r))
    r0,_=Rof(arm,lambda i,m:m['copy']==0); r1,_=Rof(arm,lambda i,m:m['copy']==1); d=r0-r1
    b=np.median(rng.choice(d,(4000,len(d))),1); print('  copy gap median %+.4f CI [%+.4f,%+.4f]  mean %+.4f  ; by-tau-window: '%(np.median(d),*np.percentile(b,[2.5,97.5]),d.mean()))
# arm 1 copy gap decomposed: per-tau median C difference
r0,C0_=Rof('1',lambda i,m:m['copy']==0); r1,C1_=Rof('1',lambda i,m:m['copy']==1)
print('arm1 copy0-copy1 median C diff by tau (26..78):',' '.join('%d:%+.3f'%(t,np.median(C0_[:,k]-C1_[:,k])) for k,t in enumerate(taus) if 26<=t<=78))
# 4: degree
Er,zz=zzmap(props)
deg={a:sum(1 for b in B if tuple(sorted((a,b))) in zz) for a in A}
for arm in ('+','1','0'):
    X,Y,C=curves(M,taus,arm); w=lambda lo,hi:C[:,(taus>=lo)&(taus<=hi)].mean(1); R=w(58,78)-w(26,50)
    for d in (2,3,4):
        m=np.array([deg[a]==d for a in A]); 
        if m.sum()>=3:
            bb=np.median(rng.choice(R[m],(4000,m.sum())),1); print('arm',arm,'deg',d,'n',m.sum(),'median R %+.4f [%+.4f,%+.4f]'%(np.median(R[m]),*np.percentile(bb,[2.5,97.5])))
X,Y,C=curves(M,taus,'+'); w=lambda lo,hi:C[:,(taus>=lo)&(taus<=hi)].mean(1); R=w(58,78)-w(26,50)
m3=np.array([deg[a]==3 for a in A]);m4=np.array([deg[a]==4 for a in A])
d=[np.median(rng.choice(R[m4],m4.sum()))-np.median(rng.choice(R[m3],m3.sum())) for _ in range(4000)]
print('deg4 - deg3 (arm +): %+.4f [%+.4f,%+.4f]'%(np.median(R[m4])-np.median(R[m3]),*np.percentile(d,[2.5,97.5])))
# mean-zz weighting
mz=np.array([np.mean([zz[tuple(sorted((a,b)))] for b in B if tuple(sorted((a,b))) in zz])*1e-3 for a in A])
from scipy.stats import spearmanr
print('arm + R vs mean edge zz (props): spearman',spearmanr(mz,R)); print('arm + R vs degree:',spearmanr([deg[a] for a in A],R))
# floor/noise
print('hypot noise floor: tail C (120-150us) arm+ median %.3f ; arm1 %.3f'%(np.median(curves(M,taus,'+')[2][:,taus>=120]),np.median(curves(M,taus,'1')[2][:,taus>=120])))
