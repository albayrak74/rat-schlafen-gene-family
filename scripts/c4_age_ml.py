import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
rng=np.random.default_rng(2026)
C=pd.read_pickle('out/bodymap_counts.pkl'); meta=pd.read_csv('out/bodymap_meta.csv',index_col=0)
ren={'Slfn1':'Slfn1','Slfn2':'Slfn2','Slfn3':'Slfn4','Slfn4':'Slfn3','Slfn5':'Slfn5','Slfn8':'Slfn9','LOC287569':'Slfn14','Slfnl1':'Slfnl1'}
L=np.log2(C/C.sum()*1e6+1).rename(index=ren)
M=meta[meta.week.isin([21,104])].copy(); M['y']=(M.week==104).astype(int)
L=L[M.index]
# genes expressed in every organ (mean log2CPM>1 within each organ)
ok=None
for o,s in M.groupby('organ'):
    e=set(L.index[(L[s.index].mean(1)>0.5)&(L[s.index].std(1)>0)]); ok=e if ok is None else ok&e
ok=sorted(ok); print('universe',len(ok))
Z=L.loc[ok].copy()
for o,s in M.groupby('organ'):
    x=L.loc[ok,s.index]; Z[s.index]=x.sub(x.mean(1),axis=0).div(x.std(1)+1e-9,axis=0)
SLFN=[g for g in ['Slfn1','Slfn2','Slfn3','Slfn4','Slfn5','Slfn9','Slfn14','Slfnl1'] if g in ok]; print('SLFN',SLFN)
ISG=[g for g in ['Irf7','Isg15','Mx1','Rsad2','Stat1','Ifit3','Oasl','Usp18'] if g in ok][:len(SLFN)]; print('ISG',ISG)
organs=sorted(M.organ.unique())
def lolo(genes,model='en'):
    ys=[];ps=[];per={}
    for o in organs:
        tr=M.index[M.organ!=o]; te=M.index[M.organ==o]
        m=LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=0.5,C=0.5,max_iter=5000) if model=='en' else RandomForestClassifier(n_estimators=300,min_samples_leaf=2,random_state=0)
        m.fit(Z.loc[genes,tr].T,M.y[tr]); p=m.predict_proba(Z.loc[genes,te].T)[:,1]
        per[o]=roc_auc_score(M.y[te],p); ys+=list(M.y[te]); ps+=list(p)
    return roc_auc_score(ys,ps),per,ps
rows=[];preds={}
for name,genes in [('Schlafen family',SLFN),('Interferon-stimulated genes',ISG)]:
    for mdl in ['en','rf']:
        a,per,ps=lolo(genes,mdl); rows.append([name,mdl,len(genes),a,*per.values()]); preds[(name,mdl)]=ps
Cn=pd.concat([L.loc[ok,s.index].sub(L.loc[ok,s.index].mean(1),axis=0) for o,s in M.groupby('organ')],axis=1)
var=Cn.var(1).sort_values(ascending=False).index[:2000]
a,per,ps=lolo(list(var),'en'); rows.append(['Transcriptome (2,000 most variable)','en',2000,a,*per.values()]); preds[('tx','en')]=ps
P=pd.DataFrame(rows,columns=['feature_set','model','k','pooled_AUC',*organs]); pd.set_option('display.width',250); print(P.round(3).to_string(index=False)); P.to_csv('out/age_ml_lolo.csv',index=False)
pd.DataFrame(preds,index=[i for o in organs for i in M.index[M.organ==o]]).to_csv('out/age_ml_preds.csv')
# bootstrap CI over animals for pooled AUC (Schlafen en)
ys=np.array([M.y[i] for o in organs for i in M.index[M.organ==o]]); ps=np.array(preds[('Schlafen family','en')])
bs=[]
for b in range(1000):
    idx=rng.integers(0,len(ys),len(ys))
    if len(set(ys[idx]))>1: bs.append(roc_auc_score(ys[idx],ps[idx]))
print('Schlafen EN pooled AUC CI',np.percentile(bs,[2.5,97.5]).round(3))
null=[]
for i in range(200):
    g=list(rng.choice(ok,len(SLFN),replace=False)); a,per,_=lolo(g,'en'); null.append([a,*per.values()])
N=pd.DataFrame(null,columns=['pooled_AUC',*organs]); N.to_csv('out/age_ml_null.csv',index=False)
obs=P.iloc[0]; print('null median',N.pooled_AUC.median().round(3),'95th',N.pooled_AUC.quantile(0.95).round(3),'P',((N.pooled_AUC>=obs.pooled_AUC).sum()+1)/(len(N)+1))
# SHAP
import shap
rf=RandomForestClassifier(n_estimators=500,min_samples_leaf=2,random_state=0).fit(Z.loc[SLFN].T,M.y)
sv=shap.TreeExplainer(rf).shap_values(Z.loc[SLFN].T); sv=sv[1] if isinstance(sv,list) else sv[:,:,1]
X=Z.loc[SLFN].T
S=pd.DataFrame({'mean_abs_shap':np.abs(sv).mean(0),'direction':[np.corrcoef(X[g],sv[:,i])[0,1] for i,g in enumerate(SLFN)]},index=SLFN); print(S.round(3)); S.to_csv('out/age_ml_shap.csv'); np.save('out/age_shap_vals.npy',sv); X.to_csv('out/age_X.csv')
# external: Wistar liver GSE160153 650d vs 110d
ex=pd.read_pickle('out/wistar_expr.pkl')['GSE160153']; lf,mt=ex
lf=lf.rename(index={'Slfn3':'Slfn4','Slfn13':'Slfn9'})
gx=[g for g in SLFN if g in lf.index]; print('external genes',gx)
zz=lf.loc[gx]; zz=zz.sub(zz.mean(1),axis=0).div(zz.std(1)+1e-9,axis=0)
m=LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=0.5,C=0.5,max_iter=5000).fit(Z.loc[gx].T,M.y)
yy=(mt.age=='650d').astype(int); print('External Wistar liver AUC',round(roc_auc_score(yy,m.predict_proba(zz.T)[:,1]),3))
