from load import *
import re
def zzmap(props):
    E=edges(props); zz={}
    for x in props['general']:
        if x['name'].startswith('zz_'):
            s=x['name'][3:]
            for k in range(1,len(s)):
                a,b=int(s[:k]),int(s[k:])
                if tuple(sorted((a,b))) in E: zz[tuple(sorted((a,b)))]=x['value']*1e9  # Hz
    return E,zz
def qp(props,name):
    r=[]
    for q in props['qubits']:
        v=[p['value'] for p in q if p['name']==name]; r.append(v[0] if v else np.nan)
    r=np.array(r,float); return np.where(np.isfinite(r),r,np.nanmedian(r))
def g1(t,z,T1):  # complex coherence factor for neighbour starting in |1>, t in s
    s=np.linspace(0,t,2001)
    integ=np.trapezoid(np.exp(-s/T1)/T1*np.exp(2j*np.pi*z*s),s,axis=-1) if False else np.trapezoid(np.exp(-s/T1)/T1*np.exp(2j*np.pi*z*s),s)
    return np.exp(-t/T1)*np.exp(2j*np.pi*z*t)+integ
def predict(chip,taus,t2mode='T2'):
    rec,props,A,B,_,M=load(chip); E,zz=zzmap(props)
    T1=qp(props,'T1')*1e-6; T2=qp(props,'T2')*1e-6
    out={'0':[], '1':[], '+':[]}; nb={}
    for a in A:
        nbs=[b for b in B if tuple(sorted((a,b))) in zz]; nb[a]=nbs
        for arm in '01+':
            c=[]
            for t in taus*1e-6:
                g=1+0j
                for b in nbs:
                    z=zz[tuple(sorted((a,b)))]; h=g1(t,z,T1[b])
                    if arm=='1': g*=h
                    elif arm=='+': g*=0.5*(1+h)
                c.append(abs(g)*np.exp(-t/T2[a]))
            out[arm].append(c)
    return A,B,nb,zz,{k:np.array(v) for k,v in out.items()}
if __name__=='__main__':
    for chip in ('phoenix','fez'):
        rec,props,A,B,taus,M=load(chip)
        A,B,nb,zz,P=predict(chip,taus)
        print(chip,'missing-edge A qubits',[a for a in A if not nb[a]][:5],'deg dist',np.bincount([len(nb[a]) for a in A]))
        for arm in '01+':
            X,Y,C=curves(M,taus,arm)
            w=lambda Z,lo,hi:Z[:,(taus>=lo)&(taus<=hi)].mean(1)
            Rm=w(C,58,78)-w(C,26,50); Rp=w(P[arm],58,78)-w(P[arm],26,50)
            print(' arm',arm,'pred median R %+.4f  meas %+.4f  corr per-qubit %.3f'%(np.median(Rp),np.median(Rm),np.corrcoef(Rp,Rm)[0,1]))
            print('   pred median C ',' '.join('%5.2f'%v for v in np.median(P[arm],0)[::3]))
            print('   meas median C ',' '.join('%5.2f'%v for v in np.median(C,0)[::3]))
