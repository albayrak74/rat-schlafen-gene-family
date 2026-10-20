import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import statsmodels.formula.api as smf
rng=np.random.default_rng(2026)
out={}
# ---------- AGE (BodyMap) ----------
CNT=pd.read_pickle('out/bodymap_counts.pkl'); meta=pd.read_csv('out/bodymap_meta.csv',index_col=0)
ren={'Slfn1':'Slfn1','Slfn2':'Slfn2','Slfn3':'Slfn4','Slfn4':'Slfn3','Slfn5':'Slfn5','Slfn8':'Slfn9','LOC287569':'Slfn14','Slfnl1':'Slfnl1'}
L=np.log2(CNT/CNT.sum()*1e6+1).rename(index=ren)
M=meta[meta.week.isin([21,104])].copy(); M['y']=(M.week==104).astype(int); L=L[M.index]
ok=None
for o,s in M.groupby('organ'):
    e=set(L.index[(L[s.index].mean(1)>0.5)&(L[s.index].std(1)>0)]); ok=e if ok is None else ok&e
ok=sorted(ok); Z=L.loc[ok].copy(); Cn=L.loc[ok].copy()
for o,s in M.groupby('organ'):
    x=L.loc[ok,s.index]; Cn[s.index]=x.sub(x.mean(1),axis=0); Z[s.index]=Cn[s.index].div(x.std(1)+1e-9,axis=0)
SL=['Slfn2','Slfn4','Slfn5','Slfn9']
def lolo(genes,MM):
    ys=[];ps=[]
    for o in sorted(MM.organ.unique()):
        tr=MM.index[MM.organ!=o]; te=MM.index[MM.organ==o]
        m=LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=0.5,C=0.5,max_iter=5000).fit(Z.loc[genes,tr].T,MM.y[tr])
        ys+=list(MM.y[te]); ps+=list(m.predict_proba(Z.loc[genes,te].T)[:,1])
    return roc_auc_score(ys,ps)
top=Cn.var(1).sort_values(ascending=False).index[:2000]
out['age_tx_withinorgan_var']=lolo(list(top),M)
Mn=M[M.organ!='Tst']
obs=lolo(SL,Mn); null=[lolo(list(rng.choice(ok,4,replace=False)),Mn) for _ in range(200)]
out['age_noTestis_SLFN']=obs; out['age_noTestis_null_median']=float(np.median(null)); out['age_noTestis_P']=(sum(n>=obs for n in null)+1)/201
null_all=pd.read_csv('out/age_ml_null.csv').pooled_AUC; out['age_null_k']=int((null_all>=0.707).sum())
# Ptprc-adjusted age effects (within organ, OLS on log2CPM)
rows=[]
for o,s in M.groupby('organ'):
    for g in ['Slfn2','Slfn4','Slfn5','Slfn9']:
        d=pd.DataFrame({'y':L.loc[g,s.index],'age':s.y,'sex':s.sex,'ptprc':L.loc['Ptprc',s.index],'adgre1':L.loc['Adgre1',s.index] if 'Adgre1' in L.index else 0})
        if d.y.mean()<0.5: continue
        f0='y~age'+('+C(sex)' if d.sex.nunique()>1 else ''); f1=f0+'+ptprc'
        a=smf.ols(f0,d).fit(); b=smf.ols(f1,d).fit()
        rows.append([o,g,a.params['age'],a.pvalues['age'],b.params['age'],b.pvalues['age']])
P=pd.DataFrame(rows,columns=['organ','gene','age_beta','age_p','age_beta_adjPtprc','age_p_adjPtprc']); P.to_csv('out/age_ptprc_adjusted.csv',index=False)
pt=[]
for o,s in M.groupby('organ'):
    a=L.loc['Ptprc',s.index]; pt.append([o,a[s.y==1].mean()-a[s.y==0].mean()])
out['ptprc_age_diff']={o:round(v,2) for o,v in pt}
# ---------- INFLAMMATION ----------
src=open('c3_ml.py').read().split('SLFN=[g')[0]
ns={}; exec(src,ns)
Zi=ns['Z']; D=ns['D']; common=ns['common']
Cv={}
for k,(m,y,bt) in D.items():
    mm=m.loc[common]
    if bt is None: Cv[k]=mm.sub(mm.mean(1),axis=0)
    else:
        o=mm.copy()
        for b_ in bt.unique(): cols=bt.index[bt==b_]; o[cols]=mm[cols].sub(mm[cols].mean(1),axis=0)
        Cv[k]=o
var=pd.concat(Cv.values(),axis=1).var(1).sort_values(ascending=False).index[:2000]
ys=[];ps=[];per={}
for test in Zi:
    Xtr=pd.concat([Zi[k][0].loc[var].T for k in Zi if k!=test]); ytr=pd.concat([Zi[k][1] for k in Zi if k!=test])
    m=LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=0.5,C=0.5,max_iter=5000).fit(Xtr.fillna(0),ytr)
    p=m.predict_proba(Zi[test][0].loc[var].T.fillna(0))[:,1]; per[test]=roc_auc_score(Zi[test][1],p); ys+=list(Zi[test][1]); ps+=list(p)
out['infl_tx_withinds_var']=roc_auc_score(ys,ps); out['infl_tx_per']={k:round(v,3) for k,v in per.items()}
import json; print(json.dumps(out,indent=1,default=float)); json.dump(out,open('out/revision_ml.json','w'),indent=1,default=float)
print(P.round(3).to_string())
