import sys;sys.path.insert(0,'/private/tmp/claude-501/-Users-heuremen/80087010-0cb3-4ae5-9e13-362e6d5086ad/scratchpad/ref_zz')
from lib import *
chip=sys.argv[1];Jb=J(chip)
# orientation sanity: b1 mean in Z arm (should be ~readout err), P1 arm (~0.5), O1 (~1)
for arm in ('Z','P1','O1'):
    i=Jb.ix[(arm,Jb.taus[3],'X',0)]; print(arm,'b1 mean per triple (first 6)',np.round(Jb.bit('b1',i).mean(0)[:6],2),'b2',np.round(Jb.bit('b2',i).mean(0)[:3],2))
taus=np.array(Jb.taus);ts=taus*1e-6
C={a:np.array([coh(Jb,a,t,(0,1)) for t in Jb.taus]) for a in ('Z','P1','P2','P12','O1','EZ','EP1')}
G=np.arange(0.10,2.5001,0.0025)
f1=np.array([Jb.zz[tuple(sorted((a,b1)))] for b1,a,b2 in Jb.tr])
s=np.full(Jb.m,np.nan)
for j in range(Jb.m):
    if f1[j]<=0:continue
    mod=C['Z'][:,j][None]*np.abs(np.cos(np.pi*G[:,None]*f1[j]*ts[None]));s[j]=G[np.argmin(((mod-C['P1'][:,j][None])**2).sum(1))]
print('R1 median',np.nanmedian(s),'n',np.isfinite(s).sum())
rng=np.random.default_rng(1);x=s[np.isfinite(s)];bs=np.median(rng.choice(x,(2000,len(x))),1);print(np.percentile(bs,[2.5,97.5]))
print('s per triple',np.round(s,2))
# R3
T1b=np.array([qp(Jb.p,b1,'T1') for b1,a,b2 in Jb.tr])
du=np.full(Jb.m,np.nan);dc=du.copy()
for j in range(Jb.m):
    if not np.isfinite(s[j]):continue
    td=0.5/(s[j]*f1[j])*1e6
    if td>taus[-2]:continue
    k=int(np.argmin(abs(taus-td)));t=Jb.taus[k]
    du[j]=C['P1'][k,j]/C['Z'][k,j]
    tot=0;w=0
    for b in(0,1):
        x_,n=expv(Jb,'P1',t,'X',(0,1),('b1',b));y_,_=expv(Jb,'P1',t,'Y',(0,1),('b1',b));tot+=np.hypot(x_,y_)[j]*n[j];w+=n[j]
    cb=tot/w
    # half-split of Z using own rng
    tot=0;w=0
    for h in(0,1):
        r2=np.random.default_rng(100+h);xx=0;yy=0;nn=0
        res=[]
        for st in 'XY':
            num=0;den=0
            for c in(0,1):
                i=Jb.ix[('Z',t,st,c)];v=1-2*Jb.bit('a',i)[:,j].astype(float);mask=r2.random(len(v))<0.5
                num+=(v*mask).sum();den+=mask.sum()
            res.append((num/den,den))
        tot+=np.hypot(res[0][0],res[1][0])*res[0][1];w+=res[0][1]
    ch=tot/w
    dc[j]=cb/ch
print('R3a',np.nanmedian(du),'R3b',np.nanmedian(dc),'R3c',np.nanmedian(dc-du),'n',np.isfinite(du).sum())
