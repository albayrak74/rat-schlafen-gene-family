import pandas as pd, numpy as np, re, warnings; warnings.filterwarnings('ignore')
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import stats
import statsmodels.formula.api as smf
G=['Slfn1','Slfn2','Slfn3','Slfn4','Slfn5','Slfn9','Slfn14','Slfnl1']
ENS={'ENSRNOG00000009906':'Slfnl1','ENSRNOG00000021412':'Slfn9','ENSRNOG00000021719':'Slfn5','ENSRNOG00000037113':'Slfn2','ENSRNOG00000057092':'Slfn4','ENSRNOG00000048053':'Slfn1'}
rows=[]; expr={}
def deseq(counts,meta,design,contrast,ds,label,genes_map):
    cc=counts.loc[:,meta.index]; cc=cc[(cc>=10).sum(1)>=3].round().astype(int)
    dds=DeseqDataSet(counts=cc.T,metadata=meta,design=design,quiet=True); dds.deseq2()
    st=DeseqStats(dds,contrast=contrast,quiet=True); st.summary(); r=st.results_df
    for gid,g in genes_map.items():
        if gid in r.index: x=r.loc[gid]; rows.append([ds,label,g,x.baseMean,x.log2FoldChange,x.lfcSE,x.pvalue,x.padj])
        else: rows.append([ds,label,g,np.nan,np.nan,np.nan,np.nan,np.nan])
    return cc
def lognorm(cc): return np.log2(cc/cc.sum()*1e6+1)
# 1 GSE149517 TNBS colitis (counts, symbols Rnor6 NCBI: Slfn3=GeneID114247->Slfn4, Slfn13->Slfn9)
x=pd.read_excel('geo/GSE149517_Gene_expressions.xlsx').set_index('Unnamed: 0'); x=x[~x.index.duplicated()]
m149={'Slfn1':'Slfn1','Slfn2':'Slfn2','Slfn3':'Slfn4','Slfn5':'Slfn5','Slfn13':'Slfn9','Slfnl1':'Slfnl1'}
meta=pd.DataFrame(index=x.columns.astype(str)); x.columns=meta.index
meta['grp']=['Inflamed' if c.endswith('gy') else ('Control' if c in ('216','217') else 'Uninflamed') for c in meta.index]
meta['animal']=[re.sub('gy','',c) for c in meta.index]
for lab,keep in [('TNBS inflamed vs control',['Control','Inflamed']),('TNBS uninflamed vs control',['Control','Uninflamed'])]:
    mm=meta[meta.grp.isin(keep)].copy(); deseq(x,mm,'~grp',['grp',keep[1],'Control'],'GSE149517',lab,m149)
mm=meta[meta.grp!='Control'].copy(); deseq(x,mm,'~animal+grp',['grp','Inflamed','Uninflamed'],'GSE149517','TNBS inflamed vs uninflamed (paired)',m149)
expr['GSE149517']=(lognorm(x),meta)
# 2 GSE228906 LPS lung FPKM (Ensembl ids)
a=pd.read_csv('geo/GSE228906_Expression.csv.gz').set_index('Gene_ID'); a=a.drop(columns='Name'); a.columns=[c.split(':')[0] for c in a.columns]
meta=pd.DataFrame(index=a.columns); meta['grp']=['Control' if c.startswith('SD') else ('LPS' if c.startswith('SM') else 'LPS+SJY') for c in meta.index]
la=np.log2(a+1)
for gid,g in ENS.items():
    if gid in la.index:
        c=la.loc[gid,meta.grp=='Control']; t=la.loc[gid,meta.grp=='LPS']; tt=stats.ttest_ind(t,c,equal_var=False)
        rows.append(['GSE228906','LPS lung injury vs control',g,a.loc[gid].mean(),t.mean()-c.mean(),np.nan,tt.pvalue,np.nan])
expr['GSE228906']=(la,meta)
# 3 GSE261056 adjuvant arthritis liver counts
b=pd.read_csv('geo/GSE261056_gene_count_matrix.txt.gz',sep='\t').set_index('gene_id')
meta=pd.DataFrame(index=b.columns); meta['grp']=[c.split('-')[0] for c in meta.index]
mm=meta[meta.grp.isin(['NC','model'])]; deseq(b,mm,'~grp',['grp','model','NC'],'GSE261056','Adjuvant arthritis liver vs control',ENS)
expr['GSE261056']=(lognorm(b),meta)
# 4 GSE298271 IL-1b microglia counts
c=pd.read_csv('geo/GSE298271_processed_raw_counts.csv.gz',encoding='utf-8-sig',low_memory=False).set_index('Id')
sc=[s for s in c.columns if s.startswith('p')]; cc=c[sc].apply(pd.to_numeric,errors='coerce').fillna(0)
meta=pd.DataFrame(index=sc); meta['grp']=['PBS' if '-pbs-' in s else ('IL1b' if re.search(r'-il1-\d+$',s) else 'IL1b+MSC') for s in sc]
meta['stage']=[s.split('-')[0] for s in sc]; meta['batch']=[s.split('-')[1] for s in sc]
mm=meta[meta.grp.isin(['PBS','IL1b'])].copy(); mm['sb']=mm.stage+'_'+mm.batch
deseq(cc,mm,'~sb+grp',['grp','IL1b','PBS'],'GSE298271','IL-1β microglia vs PBS',ENS)
expr['GSE298271']=(lognorm(cc),meta)
# 5 GSE192425 DIO obesity liver (current RefSeq names), batch A: Lean vs Saline+Sham
d=pd.read_csv('geo/GSE192425_Raw_count_table.csv.gz',sep='\t'); d=d.drop_duplicates('gene').set_index('gene')
sc=[s for s in d.columns if s.startswith('L_')]; dd=d[sc].astype(float)
meta=pd.DataFrame(index=sc); meta['grp']=[s.split('_')[1] for s in sc]; meta['batch']=[s.split('_')[-1] for s in sc]
mm=meta[(meta.batch=='batchA')&meta.grp.isin(['Lean','Saline','Sham'])].copy(); mm['ob']=np.where(mm.grp=='Lean','Lean','Obese')
deseq(dd,mm,'~ob',['ob','Obese','Lean'],'GSE192425','Diet-induced obesity liver vs lean',{**{g:g for g in G if g!='Slfn9'},'Slfn13':'Slfn9'})
expr['GSE192425']=(lognorm(dd),meta)
# 6 GSE160153 liver ageing x sex (normalised values; RefSeq 2020: Slfn3=GeneID114247->Slfn4, Slfn13->Slfn9)
fm=[]
for sx in ['MALE','FEMALE']:
    t=pd.read_csv(f'geo/GSE160153_PROCESSED_DATA_{sx}.txt.gz',sep='\t'); t=t.drop_duplicates('Title').set_index('Title'); fm.append(t[[c for c in t.columns if re.match(r'\d+d_',c)]])
f=pd.concat(fm,axis=1); lf=np.log2(f.astype(float))
meta=pd.DataFrame(index=f.columns); meta['age']=[c.split('_')[0] for c in f.columns]; meta['sex']=[c.split('_')[1] for c in f.columns]
m160={'Slfn3':'Slfn4','Slfn13':'Slfn9','Slfn2':'Slfn2','Slfnl1':'Slfnl1','Slfn1':'Slfn1','Slfn4':'Slfn3','Slfn5':'Slfn5'}
for old,g in m160.items():
    if old in lf.index:
        dfm=meta.copy(); dfm['y']=lf.loc[old].values
        fit=smf.ols('y~C(age)+C(sex)',dfm).fit()
        ka=[k for k in fit.params.index if 'age' in k][0]; ks=[k for k in fit.params.index if 'sex' in k][0]
        rows.append(['GSE160153','Liver 650 d vs 110 d (adj. sex)',g,f.loc[old].mean(),fit.params[ka],fit.bse[ka],fit.pvalues[ka],np.nan])
        rows.append(['GSE160153','Liver male vs female (adj. age)',g,f.loc[old].mean(),fit.params[ks],fit.bse[ks],fit.pvalues[ks],np.nan])
expr['GSE160153']=(lf,meta)
R=pd.DataFrame(rows,columns=['dataset','contrast','gene','baseMean','log2FC','SE','p','padj_tx'])
from statsmodels.stats.multitest import multipletests
R['fdr_family']=np.nan
for k,s in R.groupby(['dataset','contrast']):
    ok=s.p.notna(); 
    if ok.sum(): R.loc[s.index[ok],'fdr_family']=multipletests(s.p[ok],method='fdr_bh')[1]
R.to_csv('out/wistar_slfn_DE.csv',index=False); pd.to_pickle(expr,'out/wistar_expr.pkl')
pd.set_option('display.width',250)
print(R.pivot_table(index='contrast',columns='gene',values='log2FC').round(2))
print(R.pivot_table(index='contrast',columns='gene',values='fdr_family').round(3))
