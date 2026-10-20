import requests,re,pandas as pd,random
s=requests.Session(); H={'Content-Type':'application/json'}
G=pd.read_csv('out/slfn_rat_genes.csv')
motifs={'ISRE':r'GAAA[ACGT]{2}GAAA','GAS':r'TTC[ACGT]{3}GAA','NF-kB':r'GGG[AG][ACGT]T[CT][CT]CC','IRF-E':r'AA[ACGT]{2}GAAA'}
def rc(x): return x[::-1].translate(str.maketrans('ACGTN','TGCAN'))
def count(seq):
    return {m:len(re.findall(p,seq))+len(re.findall(p,rc(seq))) for m,p in motifs.items()}
rows=[]; prom={}
for _,g in G.iterrows():
    if g.strand=='+': a,b=g.start-2000,g.start-1
    else: a,b=g.end+1,g.end+2000
    seq=s.get(f'https://rest.ensembl.org/sequence/region/rattus_norvegicus/{g.chr}:{a}..{b}:1',headers={'Content-Type':'text/plain'}).text.upper()
    if g.strand=='-': seq=rc(seq)
    prom[g.gene]=seq; c=count(seq); c.update(gene=g.gene,GC=round((seq.count('G')+seq.count('C'))/len(seq)*100,1)); rows.append(c)
# background: 200 random 2-kb windows from chr10 (shuffled-position)
random.seed(1); bg=[]
for i in range(200):
    st=random.randint(5_000_000,100_000_000)
    seq=s.get(f'https://rest.ensembl.org/sequence/region/rattus_norvegicus/10:{st}..{st+1999}:1',headers={'Content-Type':'text/plain'}).text.upper()
    if seq.count('N')<100: bg.append(count(seq))
B=pd.DataFrame(bg); P=pd.DataFrame(rows).set_index('gene')
print(P); print('background mean',B.mean().round(2).to_dict(),'n',len(B))
P.to_csv('out/slfn_promoter_motifs.csv'); B.to_csv('out/promoter_background.csv',index=False)
open('out/slfn_promoters.fa','w').write(''.join(f'>{k}\n{v}\n' for k,v in prom.items()))
