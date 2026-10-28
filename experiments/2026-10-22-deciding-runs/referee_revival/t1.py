from load import *
for chip in ('phoenix','fez'):
    rec,props,A,B,taus,M=load(chip)
    print(chip,len(A),len(B),'meta0',rec['meta'][0],'taus',taus[:3],taus[-1],len(taus))
    for arm in '01+':
        res=[]
        for cp in ((0,),(1,),(0,1)):
            X,Y,C=curves(M,taus,arm,cp)
            w=lambda lo,hi:C[:,(taus>=lo)&(taus<=hi)].mean(1)
            R=w(58,78)-w(26,50); res.append(R)
        rng=np.random.default_rng(1)
        r=res[2]; b=np.median(rng.choice(r,(2000,len(r))),1)
        print(' arm',arm,'median R %+.4f [%+.4f,%+.4f]'%(np.median(r),*np.percentile(b,[2.5,97.5])),'copy0 %+.4f copy1 %+.4f'%(np.median(res[0]),np.median(res[1])))
