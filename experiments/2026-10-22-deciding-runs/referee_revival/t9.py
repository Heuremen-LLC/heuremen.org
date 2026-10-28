from model2 import *
import model, model2
rng=np.random.default_rng(9)
rec,props,A,B,taus,M=load('phoenix'); E_,zz=zzmap(props)
for sig,sc in ((0.0,0.5),(0.1,0.5),(0.2,0.5),(0.3,0.5),(0.2,0.45),(0.2,0.55)):
    Rs=[]
    for k in range(12 if sig>0 else 1):
        zz2={e:v*sc*np.exp(rng.normal(0,sig)) for e,v in zz.items()}
        model2.zzmap=lambda p,zz2=zz2:(E_,zz2); model.zzmap=model2.zzmap
        res,_,_,_=model2.run('phoenix',verbose=False); Rs.append(np.median(res['+'][0]))
    print('jitter sigma(lognormal) %.1f scale %.2f : predicted median R+ mean %+.4f  range [%+.4f,%+.4f]  (measured +0.0386)'%(sig,sc,np.mean(Rs),min(Rs),max(Rs)))
