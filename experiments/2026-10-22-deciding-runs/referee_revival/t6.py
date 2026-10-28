from model2 import *
rng=np.random.default_rng(5)
def gfac(zs,T1s,t,arm):
    g=np.ones(len(t),complex)
    for z,T1 in zip(zs,T1s):
        h=np.array([g1(tt,z,T1) for tt in t]); g*= 0.5*(1+h) if arm=='+' else h
    return abs(g)
for chip,sc in (('phoenix',0.5),('fez',1.05)):
    rec,props,A,B,taus,M=load(chip); E_,zz=zzmap(props); t=taus*1e-6; T1=qp(props,'T1')*1e-6
    X0,Y0,C0=curves(M,taus,'0'); C1=curves(M,taus,'1')[2]; Cp=curves(M,taus,'+')[2]
    f=lambda tt,amp,T,p: amp*np.exp(-(tt/T)**p)
    ext=[];nbz=[]
    for i,a in enumerate(A):
        try: po,_=curve_fit(f,t,C0[i],p0=[0.97,100e-6,1.0],bounds=([0.5,5e-6,0.5],[1.2,1e-3,2.5]))
        except Exception: po=[0.97,100e-6,1.0]
        ext.append(f(t,*po)); nbs=[b for b in B if tuple(sorted((a,b))) in zz]
        nbz.append([(zz[tuple(sorted((a,b)))]*sc,T1[b]) for b in nbs])
    ext=np.array(ext)
    def sse(assign,arm,Cm):
        s=0
        for i in range(len(A)):
            zs=nbz[assign[i]]; p=floor(ext[i]*gfac([z for z,_ in zs],[T for _,T in zs],t,arm)); s+=np.sum((p-Cm[i])**2)
        return s
    for arm,Cm in (('+',Cp),('1',C1)):
        true=sse(np.arange(len(A)),arm,Cm)
        # permutation within same degree
        deg=np.array([len(x) for x in nbz]); perms=[]
        for _ in range(60):
            asg=np.arange(len(A))
            for d in set(deg):
                idx=np.where(deg==d)[0]; asg[idx]=rng.permutation(idx)
            perms.append(sse(asg,arm,Cm))
        perms=np.array(perms)
        print(chip,'arm',arm,'SSE true edge-zz %.3f ; shuffled (same degree) mean %.3f min %.3f ; frac shuffled<=true %.3f'%(true,perms.mean(),perms.min(),(perms<=true).mean()))
