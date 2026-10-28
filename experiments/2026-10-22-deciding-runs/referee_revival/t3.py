from model2 import *
import model
rec,props,A,B,taus,M=load('phoenix')
E_,zz=zzmap(props)
for sc in (0.3,0.4,0.45,0.5,0.55,0.6):
    zz2={k:v*sc for k,v in zz.items()}
    model.zzmap=lambda p,zz2=zz2:(E_,zz2)
    import model2; model2.zzmap=model.zzmap
    res,_,_,_=model2.run('phoenix',verbose=False)
    Rp,Rm,pred,C=res['+']
    sse=np.mean((np.median(pred,0)-np.median(C,0))**2)
    print('scale',sc,'pred R+ %+.4f corr %.3f  curve mse %.5f'%(np.median(Rp),np.corrcoef(Rp,Rm)[0,1],sse), 'arm1 mse %.5f'%np.mean((np.median(res['1'][2],0)-np.median(res['1'][3],0))**2))
    if sc in (0.5,): print('  pred',' '.join('%4.2f'%v for v in np.median(pred,0)[:20])); print('  meas',' '.join('%4.2f'%v for v in np.median(C,0)[:20]))
