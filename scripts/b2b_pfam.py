import pyhmmer,gzip,pandas as pd
from Bio import SeqIO
alpha=pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile('hmm/pfam_slfn.hmm') as f: hmms=list(f)
genes=pd.read_csv('out/slfn_rat_genes.csv')
seqs=[pyhmmer.easel.TextSequence(name=g.encode(),sequence=str(next(SeqIO.parse(f'out/{g}.faa','fasta')).seq)).digitize(alpha) for g in genes.gene]
rows=[]
for top in pyhmmer.hmmsearch(hmms,seqs,bit_cutoffs='gathering'):
    q=top.query.name if isinstance(top.query.name,str) else top.query.name.decode()
    for hit in top:
        for d in hit.domains.included:
            rows.append(dict(gene=hit.name if isinstance(hit.name,str) else hit.name.decode(),domain=q,start=d.env_from,end=d.env_to,score=round(d.score,1),ievalue=d.i_evalue))
D=pd.DataFrame(rows).sort_values(['gene','start']); D.to_csv('out/slfn_rat_domains_pfam.csv',index=False); print(D.to_string(index=False))
prot=[pyhmmer.easel.TextSequence(name=r.id.encode(),description=r.description.encode(),sequence=str(r.seq).replace('*','')).digitize(alpha) for r in SeqIO.parse(gzip.open('ref/GCF_036323735.1_GRCr8_protein.faa.gz','rt'),'fasta')]
for top in pyhmmer.hmmsearch([h for h in hmms if h.name in (b'SLFN_AlbA_2','SLFN_AlbA_2'.encode())],prot,bit_cutoffs='gathering'):
    d=sorted(set((h.description if isinstance(h.description,str) else h.description.decode()).split(' isoform')[0].split(' [')[0] for h in top))
    print('proteome AlbA_2 GA hits:',len(top),d)
