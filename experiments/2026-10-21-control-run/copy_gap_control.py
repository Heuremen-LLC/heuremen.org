"""Copy gap (copy 0 minus copy 1, with interval) for every control-run reading printed as copies-disagree, using the committed analyze.py functions. Defined after the data."""
import sys,json,os,numpy as np
D=os.path.dirname(os.path.abspath(__file__))+'/'
sys.argv=['x',D]
src=open(D+'analyze.py').read().split("if __name__")[0]
exec(src)
rng=np.random.default_rng(20261029)
def gap(v0,v1,label):
    d=np.asarray(v0)-np.asarray(v1); d=d[np.isfinite(d)]
    b=np.median(rng.choice(d,(2000,len(d))),axis=1); lo,hi=np.percentile(b,[2.5,97.5])
    print('  %-52s copy0-copy1 %+.4f [%+.4f,%+.4f] n=%d -> %s'%(label,np.median(d),lo,hi,len(d),'NO SHIFT (pooled stands)' if lo<=0<=hi else 'SHIFT (CT stands)'))
for chip in ('ibm_fez','ibm_kingston','ibm_marrakesh'):
    rec=json.load(open(D+'job_%s.json'%chip)); z=np.load(D+'raw_%s.npz'%chip)
    meta=rec['meta']; S=rec['S']; P=rec['pairs']; nS=len(S); nP=len(P)
    idx=lambda **kw:[i for i,m in enumerate(meta) if all(m.get(k)==v for k,v in kw.items())]
    def expv(arm,mode,N,setting,copies,K=None):
        kw=dict(arm=arm,mode=mode,N=N,setting=setting)
        if K is not None: kw['K']=K
        a=np.concatenate([z['fin_%d'%i] for c in copies for i in idx(copy=c,**kw)])
        return np.array([(1-2*bits(a,k)).mean() for k in range(nS)])
    def phase(arm,mode,N,copies,K=None):
        x=expv(arm,mode,N,'XX',copies,K);y=expv(arm,mode,N,'YY',copies,K);return np.arctan2(y,x),np.hypot(x,y)
    print('=====',chip)
    kick={}
    for c in (0,1):
        pM,_=phase('S','M',8,[c]);pD8,_=phase('S','D',8,[c]);kick[c]=np.degrees(wrap(pM-pD8))/8
    gap(kick[0],kick[1],'R1 kick N=8')
    for K in (1,10):
        kk={}
        for c in (0,1):
            pM,_=phase('C','M',8,[c],K);pD,_=phase('C','D',8,[c],K);kk[c]=np.degrees(wrap(pM-pD))/8
        gap(kk[0],kk[1],'R2 kick K=%d'%K)
    d={c:{} for c in (0,1)}
    for c in (0,1):
        for N in (2,4,8):
            _,cM=phase('S','M',N,[c]);_,cD=phase('S','D',N,[c]);d[c][N]=cM-cD
    for N in (2,4,8): gap(d[0][N],d[1][N],'R4 coherence M-D N=%d (all shots, all 117)'%N)
    sl={c:np.array([np.polyfit([2,4,8],[d[c][N][j] for N in (2,4,8)],1)[0] for j in range(nS)]) for c in (0,1)}
    gap(sl[0],sl[1],'R4 growth slope')
    def concs(mode,N,c):
        cells=[dict() for _ in range(nP)]
        for st in ('ZZ','XX','YY','XY','YX'):
            a=np.concatenate([z['fin_%d'%i] for i in idx(arm='P',mode=mode,N=N,setting=st,copy=c)])
            for j in range(nP):
                b1=bits(a,2*j);b2=bits(a,2*j+1)
                if st=='ZZ': cells[j]['ZZ']=[np.mean((b1==x)&(b2==y)) for x,y in ((0,0),(0,1),(1,0),(1,1))]
                else: cells[j][st]=np.mean((1-2*b1)*(1-2*b2))
        return np.array([conc(cc) for cc in cells])
    for N in (2,4):
        cDp=concs('D',N,0)*0
        c0=[concs('M',N,c)-concs('D',N,c) for c in (0,1)]
        keep=(concs('D',N,0)>0)&(concs('D',N,1)>0)
        gap(c0[0][keep],c0[1][keep],'R5 concurrence M-D N=%d'%N)
