import pandas as pd, numpy as np, glob, re, warnings; warnings.filterwarnings('ignore')
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from statsmodels.stats.multitest import multipletests
fs=sorted(glob.glob('geo/bodymap/*.txt.gz'))
mats=[]
for f in fs:
    d=pd.read_csv(f,sep='\t',index_col=0); mats.append(d.iloc[:,0].rename(d.columns[0]))
C=pd.concat(mats,axis=1).fillna(0); C=C.groupby(level=0).sum().round().astype(int)
ren={'Slfn1':'Slfn1','Slfn2':'Slfn2','Slfn3':'Slfn4','Slfn4':'Slfn3','Slfn5':'Slfn5','Slfn8':'Slfn9','LOC287569':'Slfn14','Slfnl1':'Slfnl1'}
meta=pd.DataFrame([re.match(r'(\w+)_(\w)_(\d+)_(\d)',c).groups() for c in C.columns],index=C.columns,columns=['organ','sex','week','rep']); meta['week']=meta.week.astype(int)
C.to_pickle('out/bodymap_counts.pkl'); meta.to_csv('out/bodymap_meta.csv')
cpm=np.log2(C/C.sum()*1e6+1)
S=cpm.loc[list(ren)].rename(index=ren); S=S.loc[['Slfn1','Slfn2','Slfn3','Slfn4','Slfn5','Slfn9','Slfn14','Slfnl1']]
S.T.join(meta).to_csv('out/bodymap_slfn_log2cpm.csv')
print(S.T.join(meta).groupby('organ')[S.index.tolist()].mean().round(2))
# age (21 vs 104 wk) and sex (adult 21 wk) per organ
rows=[]
for org,sub in meta.groupby('organ'):
    for test in ['age','sex']:
        if test=='age': s2=sub[sub.week.isin([21,104])].copy(); s2['grp']=np.where(s2.week==104,'B','A'); design='~sex+grp' if s2.sex.nunique()>1 else '~grp'
        else:
            if sub.sex.nunique()<2: continue
            s2=sub[sub.week.isin([21,104])].copy(); s2['grp']=np.where(s2.sex=='F','B','A'); s2['age']=s2.week.astype(str); design='~age+grp'
        cc=C[s2.index]; cc=cc[(cc>=10).sum(1)>=4]
        dds=DeseqDataSet(counts=cc.T,metadata=s2,design=design,quiet=True); dds.deseq2()
        st=DeseqStats(dds,contrast=['grp','B','A'],quiet=True); st.summary(); r=st.results_df
        for old,new in ren.items():
            if old in r.index: x=r.loc[old]; rows.append([test,org,new,x.baseMean,x.log2FoldChange,x.lfcSE,x.pvalue,x.padj])
T=pd.DataFrame(rows,columns=['test','organ','gene','baseMean','log2FC','lfcSE','p','padj_tx'])
T['fdr_family']=np.nan
for t in ['age','sex']:
    m=(T.test==t)&T.p.notna(); T.loc[m,'fdr_family']=multipletests(T.loc[m,'p'],method='fdr_bh')[1]
T.to_csv('out/bodymap_slfn_age_sex.csv',index=False)
pd.set_option('display.width',220)
for t in ['age','sex']:
    print(t); print(T[T.test==t].pivot(index='organ',columns='gene',values='log2FC').round(2)); print(T[T.test==t].pivot(index='organ',columns='gene',values='fdr_family').round(3))
