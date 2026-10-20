import pandas as pd,numpy as np
from scipy import stats
P=pd.read_csv('out/slfn_promoter_motifs.csv',index_col=0); B=pd.read_csv('out/promoter_background.csv')
rng=np.random.default_rng(1); rows=[]
for m in ['ISRE','GAS','IRF-E','NF-kB']:
    obs=P[m].sum(); sims=np.array([B[m].values[rng.integers(0,len(B),8)].sum() for _ in range(100000)])
    rows.append([m,obs,round(B[m].mean()*8,2),(sims>=obs).mean(),stats.poisson.sf(obs-1,B[m].mean()*8),(P[m]>0).sum()])
R=pd.DataFrame(rows,columns=['motif','observed_8_promoters','expected','P_resampling','P_poisson','promoters_with_motif']); print(R); R.to_csv('out/promoter_enrichment.csv',index=False)
