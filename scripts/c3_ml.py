import pandas as pd, numpy as np, re, warnings; warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
rng=np.random.default_rng(2026)
ann=pd.read_csv('geo/GSE298271_processed_raw_counts.csv.gz',encoding='utf-8-sig',low_memory=False,usecols=['Id','Gene name']).dropna().drop_duplicates('Id').set_index('Id')['Gene name'].replace({'Slfn13':'Slfn9'})
ren={'Slfn13':'Slfn9','Slfn3':'Slfn4'}  # 149517 era names
def lcpm(c): return np.log2(c/c.sum()*1e6+1)
D={}
x=pd.read_excel('geo/GSE149517_Gene_expressions.xlsx').set_index('Unnamed: 0'); x=x[~x.index.duplicated()]; x.columns=x.columns.astype(str)
x=x.rename(index=ren); y=pd.Series([1 if c.endswith('gy') else 0 for c in x.columns],index=x.columns)
D['GSE149517 colon (TNBS)']=(lcpm(x),y,None)
a=pd.read_csv('geo/GSE228906_Expression.csv.gz').set_index('Gene_ID').drop(columns='Name'); a.columns=[c.split(':')[0] for c in a.columns]
a=a[[c for c in a.columns if c[:2] in ('SD','SM')]]; a=np.log2(a+1); a.index=a.index.map(ann); a=a[a.index.notna()]; a=a.groupby(level=0).mean()
D['GSE228906 lung (LPS)']=(a,pd.Series([1 if c.startswith('SM') else 0 for c in a.columns],index=a.columns),None)
b=pd.read_csv('geo/GSE261056_gene_count_matrix.txt.gz',sep='\t').set_index('gene_id'); b=b[[c for c in b.columns if c.split('-')[0] in ('NC','model')]]
b=lcpm(b); b.index=b.index.map(ann); b=b[b.index.notna()].groupby(level=0).mean()
D['GSE261056 liver (arthritis)']=(b,pd.Series([1 if c.startswith('model') else 0 for c in b.columns],index=b.columns),None)
c=pd.read_csv('geo/GSE298271_processed_raw_counts.csv.gz',encoding='utf-8-sig',low_memory=False).set_index('Id')
sc=[s for s in c.columns if s.startswith('p') and ('-pbs-' in s or re.search(r'-il1-\d+$',s))]
cc=lcpm(c[sc].apply(pd.to_numeric,errors='coerce').fillna(0)); cc.index=cc.index.map(ann); cc=cc[cc.index.notna()].groupby(level=0).mean()
D['GSE298271 microglia (IL-1β)']=(cc,pd.Series([0 if '-pbs-' in s else 1 for s in sc],index=sc),pd.Series(['-'.join(s.split('-')[:2]) for s in sc],index=sc))
# common expressed genes
common=None
for k,(m,y,bt) in D.items():
    expd=set(m.index[(m.mean(1)>0.3)&(m.std(1)>0)]); common=expd if common is None else common&expd
common=sorted(common); print('common expressed genes',len(common))
def zs(m,bt):
    if bt is None: return m.sub(m.mean(1),axis=0).div(m.std(1)+1e-9,axis=0)
    out=m.copy()
    for b_ in bt.unique():
        cols=bt.index[bt==b_]; s=m[cols]; out[cols]=s.sub(s.mean(1),axis=0).div(s.std(1)+1e-9,axis=0)
    return out
Z={k:(zs(m.loc[common],bt),y) for k,(m,y,bt) in D.items()}
for k,(z,y) in Z.items(): print(k,'n=',len(y),'pos=',int(y.sum()))
SLFN=[g for g in ['Slfn1','Slfn2','Slfn4','Slfn5','Slfn9','Slfnl1'] if g in common]; print('Slfn features',SLFN)
ISG=[g for g in ['Irf7','Isg15','Mx1','Stat1','Rsad2','Ifit3','Oasl'] if g in common][:len(SLFN)]; print('ISG',ISG)
def lodo(genes,model='en'):
    ys=[];ps=[];per={}
    for test in Z:
        Xtr=pd.concat([Z[k][0].loc[genes].T for k in Z if k!=test]); ytr=pd.concat([Z[k][1] for k in Z if k!=test])
        Xte=Z[test][0].loc[genes].T; yte=Z[test][1]
        m=LogisticRegression(penalty='elasticnet',solver='saga',l1_ratio=0.5,C=0.5,max_iter=5000) if model=='en' else RandomForestClassifier(n_estimators=300,min_samples_leaf=2,random_state=0)
        m.fit(Xtr.fillna(0),ytr); p=m.predict_proba(Xte.fillna(0))[:,1]
        per[test]=roc_auc_score(yte,p); ys+=list(yte); ps+=list(p)
    return roc_auc_score(ys,ps),per
from sklearn.model_selection import LeaveOneOut, cross_val_predict
def within(genes):
    out={}
    for k,(z,y) in Z.items():
        X=z.loc[genes].T.fillna(0); m=LogisticRegression(C=0.5,max_iter=5000)
        p=cross_val_predict(m,X,y,cv=LeaveOneOut(),method='predict_proba')[:,1]; out[k]=roc_auc_score(y,p)
    return out
W=[]
rows=[]
for name,genes in [('Schlafen family',SLFN),('Interferon-stimulated genes',ISG)]:
    for mdl in ['en','rf']:
        auc,per=lodo(genes,mdl); rows.append([name,mdl,len(genes),auc,*per.values()])
for name,genes in [('Schlafen family',SLFN),('Interferon-stimulated genes',ISG)]:
    W.append([name,*within(genes).values()])
pd.DataFrame(W,columns=['feature_set',*Z.keys()]).to_csv('out/ml_within.csv',index=False); print(pd.DataFrame(W,columns=['feature_set',*Z.keys()]).round(3).to_string(index=False))
def _cen(m,bt):
    if bt is None: return m.sub(m.mean(1),axis=0)
    o=m.copy()
    for b_ in bt.unique(): cols=bt.index[bt==b_]; o[cols]=m[cols].sub(m[cols].mean(1),axis=0)
    return o
var=pd.concat([_cen(m.loc[common],bt) for k,(m,y,bt) in D.items()],axis=1).var(1).sort_values(ascending=False).index[:2000]
auc,per=lodo(list(var),'en'); rows.append(['Transcriptome (2,000 most variable)','en',2000,auc,*per.values()])
P=pd.DataFrame(rows,columns=['feature_set','model','k','pooled_AUC',*Z.keys()]); print(P.round(3).to_string(index=False))
P.to_csv('out/ml_lodo.csv',index=False)
# random null
null=[]
for i in range(200):
    g=list(rng.choice(common,len(SLFN),replace=False)); a_,per=lodo(g,'en'); null.append([a_,*per.values()])
N=pd.DataFrame(null,columns=['pooled_AUC',*Z.keys()]); N.to_csv('out/ml_null.csv',index=False)
obs=P[(P.feature_set=='Schlafen family')&(P.model=='en')].iloc[0]
print('null median pooled',N.pooled_AUC.median().round(3),'P=',((N.pooled_AUC>=obs.pooled_AUC).sum()+1)/(len(N)+1))
for k in Z: print(k,'obs',round(obs[k],3),'null med',round(N[k].median(),3),'P',round(((N[k]>=obs[k]).sum()+1)/(len(N)+1),3))
# univariate per-gene AUC within dataset (direction-free not used: signed)
U=[]
for k,(z,y) in Z.items():
    for g in SLFN: U.append([k,g,roc_auc_score(y,z.loc[g])])
pd.DataFrame(U,columns=['dataset','gene','AUC']).to_csv('out/ml_univariate.csv',index=False)
# SHAP on RF all data
import shap
X=pd.concat([Z[k][0].loc[SLFN].T for k in Z]); Y=pd.concat([Z[k][1] for k in Z])
rf=RandomForestClassifier(n_estimators=500,min_samples_leaf=2,random_state=0).fit(X,Y)
sv=shap.TreeExplainer(rf).shap_values(X); sv=sv[1] if isinstance(sv,list) else sv[:,:,1]
S=pd.DataFrame({'mean_abs_shap':np.abs(sv).mean(0),'direction':[np.corrcoef(X[g],sv[:,i])[0,1] for i,g in enumerate(SLFN)]},index=SLFN)
S.to_csv('out/ml_shap.csv'); np.save('out/ml_shap_vals.npy',sv); X.to_csv('out/ml_X.csv'); Y.to_csv('out/ml_Y.csv')
print(S.round(3))
