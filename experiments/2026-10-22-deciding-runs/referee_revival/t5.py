from model2 import *
import model, model2
def perq(chip):
    rec,props,A,B,taus,M=load(chip); E_,zz=zzmap(props); t=taus*1e-6
    T1=qp(props,'T1')*1e-6
    X0,Y0,C0=curves(M,taus,'0')
    f=lambda tt,amp,T,p: amp*np.exp(-(tt/T)**p)
    scs=np.arange(0.2,1.81,0.05); out=[]
    cache={}
    for i,a in enumerate(A):
        try: po,_=curve_fit(f,t,C0[i],p0=[0.97,100e-6,1.0],bounds=([0.5,5e-6,0.5],[1.2,1e-3,2.5]))
        except Exception: po=[0.97,100e-6,1.0]
        ext=f(t,*po); nbs=[b for b in B if tuple(sorted((a,b))) in zz]
        Cp=curves(M,taus,'+')[2][i]; C1=curves(M,taus,'1')[2][i]
        best=None
        for s in scs:
            g=np.ones(len(t),complex); g1c=np.ones(len(t),complex)
            for b in nbs:
                z=zz[tuple(sorted((a,b)))]*s; h=np.array([g1(tt,z,T1[b]) for tt in t])
                g*=0.5*(1+h); g1c*=h
            sse=np.sum((floor(ext*abs(g))-Cp)**2)+np.sum((floor(ext*abs(g1c))-C1)**2)
            if best is None or sse<best[0]: best=(sse,s)
        out.append((a,len(nbs),best[1],np.median([zz[tuple(sorted((a,b)))] for b in nbs])*1e-3))
    return np.array(out)
for chip in ('phoenix','fez'):
    o=perq(chip); print(chip,'per-qubit best global-scale on own edges: median %.2f IQR [%.2f,%.2f] ; implied median zeta_eff (kHz, scale*props median edge) %.1f'%(np.median(o[:,2]),*np.percentile(o[:,2],[25,75]),np.median(o[:,2]*o[:,3])))
    for d in sorted(set(o[:,1])): m=o[:,1]==d; print('  degree',int(d),'n',m.sum(),'median scale %.2f'%np.median(o[m,2]))
    np.save(chip+'_perq.npy',o)
