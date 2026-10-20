import pandas as pd
G=pd.read_csv('../out/slfn_rat_genes.csv'); D=pd.read_csv('../out/slfn_rat_domains_pfam.csv')
ab={'B3R':'N','SLFN_AlbA_2':'A','SLFN_GTPase-like':'G','SLFN-g3_helicase':'H'}
grp={'Slfn1':'I','Slfn2':'I','Slfn3':'II','Slfn4':'II','Slfn5':'III','Slfn9':'III','Slfn14':'III','Slfnl1':'–'}
order=['Slfn1','Slfn2','Slfn3','Slfn4','Slfn5','Slfn9','Slfn14','Slfnl1']
rows=[]
for g in order:
    r=G[G.gene==g].iloc[0]; d=D[D.gene==g].sort_values('start')
    rows.append([g,grp[g],str(r.GeneID),f"{r.chr}:{r.start:,}–{r.end:,} ({r.strand})",str(r.exons),r.protein,str(r.aa),f"{r.MW_kDa:.1f}",f"{r.pI:.2f}",f"{r.instability:.1f}",f"{r.GRAVY:.3f}",f"{r.aliphatic:.1f}",'+'.join(ab[x] for x in d.domain)+('+(H)' if g=='Slfn14' else '')])
T2=pd.DataFrame(rows,columns=['Gene','Group','GeneID','Location (chr:start–end, strand)','Exons','Protein','Length (aa)','MW (kDa)','pI','II','GRAVY','AI','Domains'])
T2.to_csv('../out/Table2_genes.csv',index=False)
T1=pd.DataFrame([
['GSE53960','Fischer 344','11 organs','2, 6, 21 and 104 weeks of age; both sexes','320 (4 per sex, organ and age)','Read counts','Yu et al. (2014)'],
['GSE149517','Wistar','Colon','TNBS colitis: inflamed vs uninflamed (paired) and vs control','6 inflamed, 6 uninflamed (paired), 2 control','Read counts','Boros et al. (2020)'],
['GSE228906','Wistar','Lung','Intratracheal LPS vs control','5 vs 5','FPKM','–'],
['GSE261056','Wistar (male)','Liver','Adjuvant-induced arthritis vs control','4 vs 4','Read counts','–'],
['GSE298271','Wistar Han (male)','Microglia (CD11b/c+)','Neonatal IL-1β exposure vs PBS','18 vs 18','Read counts','Bokobza et al. (2026)'],
['GSE192425','Wistar','Liver','Diet-induced obese (saline and sham groups, batch 1) vs lean','10 vs 6','Read counts','Metzner et al. (2022a)'],
['GSE160153','Wistar','Liver','650 vs 110 days of age; male vs female','6 per age and sex','Normalized values','Lomas-Soria et al. (2021)'],
],columns=['Accession','Strain','Tissue','Comparison','Samples used','Data type','Original study'])
T1.to_csv('../out/Table1_datasets.csv',index=False); print(T1); print(T2)
