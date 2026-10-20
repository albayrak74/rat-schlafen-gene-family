import json,subprocess,gzip,pyhmmer,pandas as pd
from Bio import SeqIO
ref=json.load(open('out/ref_pfam.json'))
names={'PF17057':'SLFN_N','PF04326':'AlbA_2','PF21026':'SLFN_GTPase-like','PF09848':'SLFN_g3_helicase','PF13538':'UvrD_C'}
alpha=pyhmmer.easel.Alphabet.amino()
hmms=[]
for pf,nm in names.items():
    seqs=[(k,v['seq'][s-1:e]) for k,v in ref.items() for (p,s,e,_) in v['pfam'] if p==pf]
    open(f'hmm/{pf}.fa','w').write(''.join(f'>{k}\n{s}\n' for k,s in seqs))
    subprocess.run(f'mafft --quiet --localpair --maxiterate 1000 hmm/{pf}.fa > hmm/{pf}.aln',shell=True)
    with pyhmmer.easel.MSAFile(f'hmm/{pf}.aln',digital=True,alphabet=alpha) as mf: msa=mf.read()
    msa.name=nm.encode()
    b=pyhmmer.plan7.Builder(alpha); bg=pyhmmer.plan7.Background(alpha); h,_,_=b.build_msa(msa,bg); hmms.append(h)
    print(pf,nm,len(seqs),'seqs, M=',h.M)
# search rat 8 proteins
rows=[]
genes=pd.read_csv('out/slfn_rat_genes.csv')
seqs=[]
for g in genes.gene:
    r=next(SeqIO.parse(f'out/{g}.faa','fasta')); seqs.append(pyhmmer.easel.TextSequence(name=g.encode(),sequence=str(r.seq)).digitize(alpha))
for top in pyhmmer.hmmsearch(hmms,seqs,E=1e-3):
    q=top.query.name if isinstance(top.query.name,str) else top.query.name.decode()
    for hit in top:
        for d in hit.domains.included:
            rows.append(dict(gene=hit.name if isinstance(hit.name,str) else hit.name.decode(),domain=q,start=d.env_from,end=d.env_to,ievalue=d.i_evalue))
D=pd.DataFrame(rows).sort_values(['gene','start']); D.to_csv('out/slfn_rat_domains.csv',index=False); print(D.to_string(index=False))
# proteome-wide search with AlbA_2 and N-term
prot=[pyhmmer.easel.TextSequence(name=r.id.encode(),description=r.description.encode(),sequence=str(r.seq).replace('*','')).digitize(alpha) for r in SeqIO.parse(gzip.open('ref/GCF_036323735.1_GRCr8_protein.faa.gz','rt'),'fasta')]
for top in pyhmmer.hmmsearch([hmms[0],hmms[1]],prot,E=1e-5):
    descs=set()
    for hit in top: descs.add((hit.description if isinstance(hit.description,str) else hit.description.decode()).split(' [')[0])
    print('proteome hits',top.query.name, len(top), sorted(descs))
