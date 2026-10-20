import pyhmmer,gzip,pandas as pd,requests,subprocess,json
from Bio import SeqIO
alpha=pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile('hmm/pfam_slfn.hmm') as f: hmms=list(f)
prot=[pyhmmer.easel.TextSequence(name=r.id.encode(),description=r.description.encode(),sequence=str(r.seq).replace('*','')).digitize(alpha) for r in SeqIO.parse(gzip.open('ref/GCF_036323735.1_GRCr8_protein.faa.gz','rt'),'fasta')]
rows=[]
for top in pyhmmer.hmmsearch(hmms,prot,bit_cutoffs='gathering'):
    q=top.query.name if isinstance(top.query.name,str) else top.query.name.decode()
    for h in top:
        de=h.description if isinstance(h.description,str) else h.description.decode()
        rows.append(dict(pfam=q,protein=h.name if isinstance(h.name,str) else h.name.decode(),description=de,score=round(h.score,1)))
P=pd.DataFrame(rows); P['product']=P.description.str.replace(r' isoform.*','',regex=True).str.replace(r' \[.*','',regex=True)
S=P.groupby(['pfam','product']).agg(n_proteins=('protein','nunique'),max_score=('score','max')).reset_index()
S.to_csv('out/proteome_pfam_hits.csv',index=False); print(S.to_string(index=False))
# tblastn mouse+human SLFN vs rat chr10 cluster region
seq=requests.get('https://rest.ensembl.org/sequence/region/rattus_norvegicus/10:68200000..68800000:1',headers={'Content-Type':'text/plain'}).text
open('rev/region.fa','w').write('>chr10_68200000\n'+seq+'\n')
ref=json.load(open('out/ref_pfam.json'))
open('rev/q.fa','w').write(''.join(f'>{k}\n{v["seq"]}\n' for k,v in ref.items() if 'SLFNL1' not in k.upper()))
subprocess.run('tblastn -query rev/q.fa -subject rev/region.fa -evalue 1e-5 -outfmt "6 qseqid sstart send evalue bitscore pident length" > rev/tbn.tsv',shell=True)
T=pd.read_csv('rev/tbn.tsv',sep='\t',header=None,names=['q','s','e','ev','bits','pid','len'])
T['pos']=68200000+T[['s','e']].min(1)-1
G=pd.read_csv('out/slfn_rat_genes.csv'); G=G[G.chr.astype(str)=='10']
def which(p):
    for _,g in G.iterrows():
        if g.start-2000<=p<=g.end+2000: return g.gene
    return 'OUTSIDE'
T['gene']=T.pos.map(which); print(T.groupby('gene').agg(n=('q','size'),best=('ev','min')))
print(T[T.gene=='OUTSIDE'].sort_values('ev').head(20).to_string())
T.to_csv('out/tblastn_cluster.csv',index=False)
