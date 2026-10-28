from load import *
for chip in ('phoenix','fez'):
    rec,props,A,B,taus,M=load(chip)
    print(chip)
    print('tau   '+' '.join('%5d'%t for t in taus))
    for arm in '01+':
        X,Y,C=curves(M,taus,arm); print('C',arm,' '.join('%5.2f'%v for v in np.median(C,0)))
    for arm in '01+':
        X,Y,C=curves(M,taus,arm); print('Xm',arm,' '.join('%5.2f'%v for v in np.mean(X,0)))
    for arm in '01+':
        X,Y,C=curves(M,taus,arm); print('Ym',arm,' '.join('%5.2f'%v for v in np.mean(Y,0)))
