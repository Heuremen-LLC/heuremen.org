import model, model2
from model2 import *
for chip,scs in (('fez',(0.25,0.5,0.75,1.0,1.25,1.5,2.0)),('phoenix',(0.4,0.45,0.5,0.55))):
    rec,props,A,B,taus,M=load(chip); E_,zz=zzmap(props)
    for sc in scs:
        zz2={k:v*sc for k,v in zz.items()}
        model2.zzmap=lambda p,zz2=zz2:(E_,zz2); model.zzmap=model2.zzmap
        res,_,_,_=model2.run(chip,verbose=False)
        print(chip,sc,' '.join('arm%s mse %.5f corrR %.2f'%(a,np.mean((np.median(res[a][2],0)-np.median(res[a][3],0))**2),np.corrcoef(res[a][0],res[a][1])[0,1]) for a in '1+'))
