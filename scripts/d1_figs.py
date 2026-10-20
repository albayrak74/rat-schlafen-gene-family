import pandas as pd, numpy as np, json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon
from Bio import Phylo
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.edgecolor':'#c3c2b7','xtick.color':'#52514e','ytick.color':'#52514e','axes.spines.top':False,'axes.spines.right':False})
DIV=LinearSegmentedColormap.from_list('d',['#104281','#5598e7','#f0efec','#e66767','#a31f1f'])
SEQ=LinearSegmentedColormap.from_list('s',['#f7f6f2','#a9cbf2','#2a78d6','#0d3a75'])
BLUE,ORANGE,AQUA,GREY,PURP='#2a78d6','#eb6834','#1baf7a','#898781','#7a5bd1'
GENES=['Slfn1','Slfn2','Slfn3','Slfn4','Slfn5','Slfn9','Slfn14','Slfnl1']
ORG={'Adr':'Adrenal','Brn':'Brain','Hrt':'Heart','Kdn':'Kidney','Lng':'Lung','Lvr':'Liver','Msc':'Muscle','Spl':'Spleen','Thm':'Thymus','Tst':'Testis','Utr':'Uterus'}
def save(f,n): f.savefig(f'fig/{n}.tif',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'}); f.savefig(f'fig/{n}.png',dpi=200,bbox_inches='tight'); plt.close(f)
def stars(q): return '' if pd.isna(q) else '***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else ''
DOMC={'B3R':'#1baf7a','SLFN_AlbA_2':'#2a78d6','SLFN_GTPase-like':'#eb6834','SLFN-g3_helicase':'#7a5bd1'}
DOML={'B3R':'SLFN N-terminal (B3R/PF17057)','SLFN_AlbA_2':'AlbA_2 core (PF04326)','SLFN_GTPase-like':'GTPase-like (PF21026)','SLFN-g3_helicase':'Group 3 helicase (PF09848)'}
G=pd.read_csv('out/slfn_rat_genes.csv').set_index('gene')
# ---------- Fig1 genomic organisation & synteny ----------
S=pd.read_csv('out/slfn_synteny_coords.csv',dtype={'chr':str})
f,ax=plt.subplots(figsize=(7.2,3.6)); ax.axis('off')
tracks=[('Rat chr 10 (GRCr8)',[(g,G.loc[g,'start'],G.loc[g,'end'],1 if G.loc[g,'strand']=='+' else -1) for g in GENES if g!='Slfnl1'],2.0,'Rno'),
        ('Mouse chr 11 (GRCm39)',[(r.gene,r.start,r.end,r.strand) for r in S[(S.species=='mus_musculus')&(S.chr=='11')].itertuples()],1.0,'Mmu'),
        ('Human chr 17 (GRCh38)',[(r.gene,r.start,r.end,r.strand) for r in S[(S.species=='homo_sapiens')&(S.chr=='17')].itertuples()],0.0,'Hsa')]
col={'Slfn5':'#2a78d6','SLFN5':'#2a78d6','Slfn9':'#7a5bd1','Slfn8':'#7a5bd1','SLFN11':'#7a5bd1','SLFN13':'#7a5bd1','Slfn14':'#a31f1f','SLFN14':'#a31f1f','Slfn1':'#1baf7a','Slfn2':'#1baf7a','SLFN12':'#1baf7a','SLFN12L':'#1baf7a','Slfn3':'#eb6834','Slfn4':'#eb6834'}
pos={}
for lab,gs,y,sp in tracks:
    s0=min(g[1] for g in gs); e0=max(g[2] for g in gs); span=e0-s0
    ax.plot([0.05,0.95],[y,y],color='#c3c2b7',lw=1.5,zorder=1); ax.text(0.0,y+0.36,lab,fontsize=7.5,fontweight='bold',color='#0b0b0b')
    for g,s,e,st in sorted(gs,key=lambda x:x[1]):
        x0=0.05+0.9*(s-s0)/span; x1=0.05+0.9*(e-s0)/span; w=max(x1-x0,0.012)
        c=col.get(g,GREY)
        if st==1: poly=[(x0,y-0.1),(x0+w*0.75,y-0.1),(x0+w,y),(x0+w*0.75,y+0.1),(x0,y+0.1)]
        else: poly=[(x0+w,y-0.1),(x0+w*0.25,y-0.1),(x0,y),(x0+w*0.25,y+0.1),(x0+w,y+0.1)]
        ax.add_patch(Polygon(poly,closed=True,color=c,zorder=3))
        ax.text(x0+w/2,(y-0.16-(0.13 if g in ('SLFN13',) else 0)) if sp!='Rno' else y+0.13,g,ha='center',va='top' if sp!='Rno' else 'bottom',fontsize=6.3,fontstyle='italic',rotation=0)
        pos[(sp,g)]=(x0+w/2,y)
    ax.text(0.27 if sp=='Rno' else 0.32,y+0.36,f'(cluster span {span/1e3:.0f} kb)',fontsize=6.2,color='#52514e')
for a,b in [('Slfn5','Slfn5'),('Slfn9','Slfn9'),('Slfn9','Slfn8'),('Slfn2','Slfn2'),('Slfn1','Slfn1'),('Slfn4','Slfn4'),('Slfn3','Slfn3'),('Slfn14','Slfn14')]:
    if ('Rno',a) in pos and ('Mmu',b) in pos:
        (x1,y1),(x2,y2)=pos[('Rno',a)],pos[('Mmu',b)]; ax.plot([x1,x2],[y1-0.1,y2+0.1],color='#c3c2b7',lw=0.7,zorder=0)
ax.text(0.0,-0.68,'Slfnl1 lies outside the cluster in each species (rat chr 5, mouse chr 4, human chr 1).',fontsize=6.8,color='#52514e')
ax.set_xlim(-0.01,1.0); ax.set_ylim(-1.3,2.5)
for i,(t,c) in enumerate([('Group I (Slfn1/2, SLFN12/12L)','#1baf7a'),('Group II (Slfn3/4)','#eb6834'),('Group III (Slfn5)','#2a78d6'),('Group III (Slfn8/9, SLFN11/13)','#7a5bd1'),('Group III (Slfn14)','#a31f1f')]):
    ax.add_patch(Rectangle((0.0+(i%3)*0.33,-1.0-(i//3)*0.17),0.015,0.09,color=c)); ax.text(0.02+(i%3)*0.33,-0.955-(i//3)*0.17,t,fontsize=6.2,va='center')
save(f,'Fig1')
# ---------- Fig2 phylogeny + domain architecture ----------
tr=Phylo.read('tree/slfn.contree','newick'); tr.root_with_outgroup('Hsa_SLFNL1'); tr.ladderize()
D=pd.read_csv('out/slfn_rat_domains_pfam.csv')
f=plt.figure(figsize=(7.3,5.4)); ax=f.add_axes([0.02,0.06,0.5,0.9])
terms=tr.get_terminals()
def lab(cl):
    if cl.is_terminal():
        sp,g=cl.name.split('_',1); return {'Rno':'Rat','Mmu':'Mouse','Hsa':'Human'}[sp]+' '+g
    try: return '' if cl.confidence is None or cl.confidence<70 else str(int(cl.confidence))
    except: return ''
Phylo.draw(tr,axes=ax,do_show=False,label_func=lambda c: lab(c) if c.is_terminal() else '',show_confidence=False)
for cl in tr.get_nonterminals():
    if cl.confidence and cl.confidence>=70 and cl!=tr.root:
        x=tr.distance(cl); ys=[ax.get_yticks()]
for t in ax.texts:
    s=t.get_text().strip()
    if s.startswith('Rat'): t.set_color('#a31f1f'); t.set_fontweight('bold')
    t.set_fontsize(6.5)
# bootstrap labels
yy={}
def ypos(cl):
    if cl.is_terminal(): return terms.index(cl)+1
    return np.mean([ypos(c) for c in cl.clades])
for cl in tr.get_nonterminals():
    if cl.confidence is not None and cl.confidence>=70 and cl.clades and cl!=tr.root:
        ax.text(tr.distance(cl)-0.02,ypos(cl)-0.25,str(int(cl.confidence)),fontsize=5.2,color='#52514e',ha='right')
ax.axis('off'); ax.set_title('a',loc='left',fontweight='bold',fontsize=9)
ax.plot([0.1,0.6],[len(terms)+0.8]*2,color='k',lw=0.8,clip_on=False); ax.text(0.35,len(terms)+0.6,'0.5 substitutions/site',ha='center',fontsize=6); ax.set_ylim(len(terms)+1.2,0.2)
ax2=f.add_axes([0.58,0.06,0.4,0.9]); ax2.set_xlim(0,1010); ax2.set_ylim(len(GENES)+0.5,0.3); ax2.axis('off'); ax2.set_title('b  Rat Schlafen proteins',loc='left',fontweight='bold',fontsize=9)
order=['Slfn5','Slfn9','Slfn14','Slfn1','Slfn2','Slfn3','Slfn4','Slfnl1']
for i,g in enumerate(order,1):
    L_=G.loc[g,'aa']; ax2.add_patch(Rectangle((0,i-0.06),L_,0.12,color='#c3c2b7'))
    for r in D[D.gene==g].itertuples(): ax2.add_patch(Rectangle((r.start,i-0.2),r.end-r.start,0.4,color=DOMC[r.domain]))
    if g=='Slfn14': ax2.add_patch(Rectangle((578,i-0.2),707-578,0.4,facecolor='white',edgecolor=DOMC['SLFN-g3_helicase'],hatch='////',lw=0.8))
    ax2.text(-15,i,g,ha='right',va='center',fontstyle='italic',fontsize=7); ax2.text(L_+12,i,f'{L_} aa',va='center',fontsize=6,color='#52514e')
for j,(k,c) in enumerate(DOMC.items()):
    ax2.add_patch(Rectangle((0,len(GENES)+0.15+j*0.0,),1,1,alpha=0))
    ax2.add_patch(Rectangle((20+(j%2)*470,len(GENES)-0.3+0.0+(j//2)*0.35+0.6),40,0.22,color=c)); ax2.text(70+(j%2)*470,len(GENES)+0.41+(j//2)*0.35,DOML[k],fontsize=5.8,va='center')
ax2.add_patch(Rectangle((20,len(GENES)+0.95),40,0.22,facecolor='white',edgecolor=DOMC['SLFN-g3_helicase'],hatch='////',lw=0.8)); ax2.text(70,len(GENES)+1.06,'Helicase-like (family profile only)',fontsize=5.8,va='center')
ax2.set_ylim(len(GENES)+1.4,0.3)
ax2.plot([700,900],[len(GENES)+1.25]*2,color='k',lw=0.8); ax2.text(800,len(GENES)+1.15,'200 aa',ha='center',fontsize=6)
save(f,'Fig2')
# ---------- Fig3 gene structure + promoter motifs + Ka/Ks ----------
EX=json.load(open('out/slfn_exons.json')); P=pd.read_csv('out/slfn_promoter_motifs.csv',index_col=0); B=pd.read_csv('out/promoter_background.csv'); K=pd.read_csv('out/slfn_kaks.csv')
f=plt.figure(figsize=(7.3,6.0)); gs=f.add_gridspec(2,2,height_ratios=[1.15,1],hspace=0.45,wspace=0.35)
ax=f.add_subplot(gs[0,:]); ax.set_title('a  Exon–intron structure (5′→3′)',loc='left',fontweight='bold',fontsize=9)
for i,g in enumerate(order):
    e=EX[g]; st=EX[g]['strand']; ex=e['exons']; cd=e['cds']; s0=min(a for a,b in ex); e0=max(b for a,b in ex)
    def tx(p): return (p-s0) if st=='+' else (e0-p)
    ax.plot([0,(e0-s0)/1000],[i,i],color='#898781',lw=0.8)
    for a,b in ex: x0,x1=sorted([tx(a),tx(b)]); ax.add_patch(Rectangle((x0/1000,i-0.12),(x1-x0)/1000,0.24,color='#c3c2b7'))
    for a,b in cd: x0,x1=sorted([tx(a),tx(b)]); ax.add_patch(Rectangle((x0/1000,i-0.2),(x1-x0)/1000,0.4,color=BLUE))
    ax.text(-0.6,i,g,ha='right',va='center',fontstyle='italic',fontsize=7); ax.text((e0-s0)/1000+0.4,i,f'{len(ex)} exons',va='center',fontsize=6,color='#52514e')
ax.set_ylim(len(order)-0.4,-0.6); ax.set_yticks([]); ax.set_xlabel('Gene length (kb)'); ax.spines['left'].set_visible(False)
ax.add_patch(Rectangle((40,7.1),1.2,0.3,color=BLUE)); ax.text(41.6,7.25,'CDS',fontsize=6,va='center'); ax.add_patch(Rectangle((44,7.1),1.2,0.3,color='#c3c2b7')); ax.text(45.6,7.25,'UTR',fontsize=6,va='center')
ax=f.add_subplot(gs[1,0]); ax.set_title('b  Promoter motifs (2 kb upstream)',loc='left',fontweight='bold',fontsize=9)
mot=['ISRE','GAS','IRF-E','NF-kB']; Pm=P.loc[order,mot]
im=ax.imshow(Pm.values,cmap=SEQ,aspect='auto',vmin=0,vmax=5)
for i in range(Pm.shape[0]):
    for j in range(Pm.shape[1]): ax.text(j,i,int(Pm.iat[i,j]),ha='center',va='center',fontsize=6.5,color='white' if Pm.iat[i,j]>=3 else '#0b0b0b')
ax.set_xticks(range(4),[m.replace('NF-kB','NF-κB') for m in mot],fontsize=7); ax.set_yticks(range(len(order)),order,fontstyle='italic'); ax.tick_params(length=0); ax.spines[:].set_visible(False)
ax.set_xlabel('Background mean: '+', '.join(f"{m.replace('NF-kB','NF-κB')} {__import__('decimal').Decimal(str(round(B[m].mean(),4))).quantize(__import__('decimal').Decimal('0.01'),rounding='ROUND_HALF_UP')}" for m in mot),fontsize=5.8,color='#52514e')
ax=f.add_subplot(gs[1,1]); ax.set_title('c  Selection (Ka/Ks)',loc='left',fontweight='bold',fontsize=9)
Ko=K[(K.type=='orthologue')&K.Ka_Ks.notna()].copy(); Ko['lab']=Ko.gene1.str.replace('Rno_','')+' / '+Ko.gene2.str.replace('Mmu_','m')
Kp=K[(K.type=='paralogue')&K.Ka_Ks.notna()&(K.Ks<1.0)].copy(); Kp['lab']=Kp.gene1.str.replace('Rno_','')+' / '+Kp.gene2.str.replace('Rno_','')
KK=pd.concat([Ko,Kp]); yv=np.arange(len(KK))
ax.barh(yv,KK.Ka_Ks,color=[BLUE if t=='orthologue' else ORANGE for t in KK.type],height=0.65)
ax.set_yticks(yv,KK.lab,fontsize=6.3,fontstyle='italic'); ax.invert_yaxis(); ax.axvline(1,color='#a31f1f',lw=0.8,ls='--'); ax.set_xlabel('Ka/Ks (dashed line: neutral expectation)'); ax.set_xlim(0,1.1)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE,label='rat–mouse pairs at corresponding positions'),Patch(color=ORANGE,label='rat paralogues')],frameon=False,fontsize=6,loc='upper center',bbox_to_anchor=(0.4,-0.2),ncol=2)
save(f,'Fig3')
# ---------- Fig4 BodyMap expression atlas ----------
L=pd.read_csv('out/bodymap_slfn_log2cpm.csv',index_col=0); tau=pd.read_csv('out/bodymap_tau.csv',index_col=0).iloc[:,0]
M=L.groupby('organ')[GENES].mean().rename(index=ORG)
f=plt.figure(figsize=(7.3,6.4)); gs=f.add_gridspec(2,3,height_ratios=[1,1],hspace=0.55,wspace=0.8)
ax=f.add_subplot(gs[0,:2]); im=ax.imshow(M.T.values,cmap=SEQ,aspect='auto',vmin=0,vmax=8.5)
for i in range(M.shape[1]):
    for j in range(M.shape[0]): v=M.T.iat[i,j]; ax.text(j,i,f'{v:.1f}',ha='center',va='center',fontsize=5.6,color='white' if v>5.5 else '#0b0b0b')
ax.set_yticks(range(8),[f'$\\it{{{g}}}$  (τ={tau[g]:.2f})' for g in GENES],fontsize=6.5); ax.set_xticks(range(len(M)),M.index,rotation=40,ha='right'); ax.tick_params(length=0); ax.spines[:].set_visible(False)
cb=f.colorbar(im,ax=ax,fraction=0.03,pad=0.02); cb.set_label('log$_2$(CPM+1)',fontsize=6.5); cb.outline.set_visible(False)
ax.set_title('a  Mean expression across organs (F344 rat BodyMap)',loc='left',fontweight='bold',fontsize=9)
ax=f.add_subplot(gs[0,2]); sub=L.groupby('week')[GENES].mean()
for g,c in zip(['Slfn4','Slfn9','Slfn2','Slfn5'],[ORANGE,PURP,AQUA,BLUE]): ax.plot(range(4),sub[g],marker='o',ms=3.5,color=c,lw=1.6,label=g)
ax.set_xticks(range(4),['2','6','21','104']); ax.set_xlabel('Age (weeks)'); ax.set_ylabel('Mean log$_2$(CPM+1), all organs'); ax.set_ylim(1.2,7.6); ax.legend(frameon=False,ncol=2,loc='upper left',prop={'style':'italic','size':6}); ax.set_title('b',loc='left',fontweight='bold',fontsize=9)
for k,(g,c) in enumerate([('Slfn4',ORANGE),('Slfn9',PURP),('Slfnl1',GREY)]):
    ax=f.add_subplot(gs[1,k]); t=L.groupby(['organ','week'])[g].mean().unstack()
    for o in t.index:
        hl=(o in ('Tst','Lng','Kdn','Spl')) if g!='Slfnl1' else (o=='Tst')
        ax.plot(range(4),t.loc[o],color=(c if g!='Slfnl1' else '#0b0b0b') if hl else '#c3c2b7',lw=1.4 if hl else 0.8,marker='o' if hl else None,ms=2.5)
        if hl: ax.text(3.08,t.loc[o].iloc[-1]+(0.12 if o=='Tst' and g=='Slfn9' else (-0.12 if o=='Kdn' and g=='Slfn9' else 0)),ORG[o],fontsize=5.6,va='center',color=c if g!='Slfnl1' else '#0b0b0b')
    ax.set_xticks(range(4),['2','6','21','104']); ax.set_xlabel('Age (weeks)'); ax.set_xlim(-0.2,3.9); ax.set_title(f'{"cde"[k]}  {g}',loc='left',fontweight='bold',fontsize=9,fontstyle='italic')
    if k==0: ax.set_ylabel('log$_2$(CPM+1)')
save(f,'Fig4')
# ---------- Fig5 age and sex effects across organs ----------
T=pd.read_csv('out/bodymap_slfn_age_sex.csv'); gg=['Slfn1','Slfn2','Slfn4','Slfn5','Slfn9','Slfnl1']
f,axs=plt.subplots(1,2,figsize=(7.3,3.6),gridspec_kw={'width_ratios':[1,0.85]})
for ax,t,ttl in [(axs[0],'age','a  Aged (104 wk) vs adult (21 wk)'),(axs[1],'sex','b  Female vs male')]:
    X=T[T.test==t].pivot(index='organ',columns='gene',values='log2FC').reindex(columns=gg); Q=T[T.test==t].pivot(index='organ',columns='gene',values='fdr_family').reindex(columns=gg)
    X.index=[ORG[o] for o in X.index]
    im=ax.imshow(X.values,cmap=DIV,norm=TwoSlopeNorm(0,-3,3),aspect='auto')
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            v=X.iat[i,j]
            if pd.isna(v): ax.add_patch(Rectangle((j-0.5,i-0.5),1,1,color='#e1e0d9')); ax.text(j,i,'n.e.',ha='center',va='center',fontsize=5,color='#898781')
            else: ax.text(j,i,f'{v:.2f}{stars(Q.iat[i,j])}',ha='center',va='center',fontsize=5.5,color='white' if abs(v)>2 else '#0b0b0b')
    ax.set_xticks(range(len(gg)),gg,fontstyle='italic',rotation=30); ax.set_yticks(range(len(X)),X.index); ax.tick_params(length=0); ax.spines[:].set_visible(False); ax.set_title(ttl,loc='left',fontweight='bold',fontsize=9)
cb=f.colorbar(im,ax=axs,fraction=0.025,pad=0.02); cb.set_label('log$_2$ fold change',fontsize=6.5); cb.outline.set_visible(False)
save(f,'Fig5')
# ---------- Fig6 Wistar models ----------
R=pd.read_csv('out/wistar_slfn_DE.csv')
order_c=['TNBS inflamed vs control','TNBS inflamed vs uninflamed (paired)','TNBS uninflamed vs control','LPS lung injury vs control','Adjuvant arthritis liver vs control','IL-1β microglia vs PBS','Diet-induced obesity liver vs lean','Liver 650 d vs 110 d (adj. sex)','Liver male vs female (adj. age)']
labs=['Colitis (TNBS), inflamed vs control','Colitis (TNBS), inflamed vs uninflamed','Colitis (TNBS), uninflamed vs control','Acute lung injury (LPS) vs control','Adjuvant arthritis, liver vs control','IL-1β-exposed microglia vs PBS','Diet-induced obesity, liver vs lean','Liver, 650 d vs 110 d','Liver, male vs female']
X=R.pivot_table(index='contrast',columns='gene',values='log2FC').reindex(index=order_c,columns=GENES); Q=R.pivot_table(index='contrast',columns='gene',values='fdr_family').reindex(index=order_c,columns=GENES)
f,ax=plt.subplots(figsize=(7.0,3.9)); im=ax.imshow(X.values,cmap=DIV,norm=TwoSlopeNorm(0,-4.5,4.5),aspect='auto')
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        v=X.iat[i,j]
        if pd.isna(v): ax.add_patch(Rectangle((j-0.5,i-0.5),1,1,color='#e1e0d9')); ax.text(j,i,'n.d.',ha='center',va='center',fontsize=5.5,color='#898781')
        else: ax.text(j,i,f'{v:.2f}{stars(Q.iat[i,j])}',ha='center',va='center',fontsize=6,color='white' if abs(v)>2.5 else '#0b0b0b')
ax.set_xticks(range(8),GENES,fontstyle='italic'); ax.set_yticks(range(len(labs)),labs); ax.tick_params(length=0); ax.spines[:].set_visible(False)
for y in [2.5,5.5,6.5]: ax.axhline(y,color='white',lw=2)
cb=f.colorbar(im,ax=ax,fraction=0.03,pad=0.02); cb.set_label('log$_2$ fold change',fontsize=6.5); cb.outline.set_visible(False)
save(f,'Fig6')
# ---------- Fig7 machine learning ----------
A=pd.read_csv('out/age_ml_lolo.csv'); N=pd.read_csv('out/age_ml_null.csv'); Sh=pd.read_csv('out/age_ml_shap.csv',index_col=0); Lo=pd.read_csv('out/ml_lodo.csv'); Nl=pd.read_csv('out/ml_null.csv')
f=plt.figure(figsize=(7.3,5.8)); gs=f.add_gridspec(2,2,hspace=0.55,wspace=0.4)
ax=f.add_subplot(gs[0,0]); orgs=list(ORG)
for name,c,mk in [('Schlafen family',ORANGE,'o'),('Interferon-stimulated genes',BLUE,'s'),('Transcriptome (2,000 most variable)','#0b0b0b','^')]:
    r=A[(A.feature_set==name)&(A.model=='en')].iloc[0]; ax.plot(range(len(orgs)),[r[o] for o in orgs],marker=mk,ms=4,lw=0,color=c,label=f"{name.split(' (')[0]}: {r.pooled_AUC:.2f}")
ax.axhline(0.5,color=GREY,lw=0.7,ls='--'); ax.set_xticks(range(len(orgs)),[ORG[o] for o in orgs],rotation=45,ha='right',fontsize=6.5); ax.set_ylabel('AUC (held-out organ)'); ax.set_ylim(0,1.05)
ax.legend(frameon=False,fontsize=5.8,loc='lower left',title='Pooled AUC',title_fontsize=6); ax.set_title('a  Aged vs adult, leave-one-organ-out',loc='left',fontweight='bold',fontsize=8.5)
ax=f.add_subplot(gs[0,1]); ax.hist(N.pooled_AUC,bins=25,color='#c3c2b7'); o=A[(A.feature_set=='Schlafen family')&(A.model=='en')].pooled_AUC.iloc[0]
ax.axvline(o,color=ORANGE,lw=2); ax.text(o+0.005,ax.get_ylim()[1]*0.9,f'Schlafen\n{o:.2f}',color=ORANGE,fontsize=6.5)
ax.set_xlabel('Pooled AUC of random 4-gene sets'); ax.set_ylabel('Count'); ax.set_title('b  Comparison with 200 random gene sets',loc='left',fontweight='bold',fontsize=8.5)
ax=f.add_subplot(gs[1,0]); sh=Sh.sort_values('mean_abs_shap'); ax.barh(range(len(sh)),sh.mean_abs_shap,color=[ORANGE if d>0.1 else (BLUE if d<-0.1 else GREY) for d in sh.direction])
ax.set_yticks(range(len(sh)),sh.index,fontstyle='italic'); ax.set_xlabel('Mean |SHAP| (random forest)'); ax.set_title('c  Contribution to the age model',loc='left',fontweight='bold',fontsize=8.5)
ax.legend(handles=[Patch(color=ORANGE,label='higher expression → aged'),Patch(color=BLUE,label='higher expression → adult'),Patch(color=GREY,label='no consistent direction')],frameon=False,fontsize=5.8,loc='lower right')
ax=f.add_subplot(gs[1,1]); dsn=[c for c in Lo.columns if c.startswith('GSE')]
for name,c,mk in [('Schlafen family',ORANGE,'o'),('Interferon-stimulated genes',BLUE,'s'),('Transcriptome (2,000 most variable)','#0b0b0b','^')]:
    r=Lo[(Lo.feature_set==name)&(Lo.model=='en')].iloc[0]; ax.plot(range(len(dsn)),[r[d] for d in dsn],marker=mk,ms=4.5,lw=0,color=c,label=f"{name.split(' (')[0]}: {r.pooled_AUC:.2f}")
ax.axhline(0.5,color=GREY,lw=0.7,ls='--'); ax.set_xticks(range(len(dsn)),['Colitis\n(colon)','LPS\n(lung)','Arthritis\n(liver)','IL-1β\n(microglia)'],fontsize=6.5); ax.set_ylabel('AUC (held-out model)'); ax.set_ylim(-0.03,1.38); ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
ax.legend(frameon=False,fontsize=5.8,loc='upper left',title='Pooled AUC',title_fontsize=6,ncol=1,borderaxespad=0.1); ax.set_title('d  Inflamed vs control, leave-one-dataset-out',loc='left',fontweight='bold',fontsize=8.5)
save(f,'Fig7')
print('ok')
