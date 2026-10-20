import gzip,re,pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
G={}; M={}; CDS={}; EX={}
chrname={}
with gzip.open('ref/GCF_036323735.1_GRCr8_genomic.gff.gz','rt') as f:
    for l in f:
        if l.startswith('#'): continue
        p=l.rstrip('\n').split('\t')
        if len(p)<9: continue
        a=dict(x.split('=',1) for x in p[8].split(';') if '=' in x)
        if p[2]=='region' and 'chromosome' in a: chrname[p[0]]=a['chromosome']
        if p[2]=='gene' and a.get('Name','').startswith('Slfn'):
            G[a['ID']]=dict(gene=a['Name'],seqid=p[0],start=int(p[3]),end=int(p[4]),strand=p[6],geneid=re.search(r'GeneID:(\d+)',a['Dbxref']).group(1))
        elif p[2]=='mRNA' and a.get('Parent') in G:
            M[a['ID']]=a['Parent']
        elif p[2]=='exon' and a.get('Parent') in M:
            EX.setdefault(a['Parent'],[]).append((int(p[3]),int(p[4])))
        elif p[2]=='CDS' and a.get('Parent') in M:
            CDS.setdefault(a['Parent'],{'prot':a.get('protein_id'),'len':0}); CDS[a['Parent']]['len']+=int(p[4])-int(p[3])+1
prot={r.id:str(r.seq) for r in SeqIO.parse(gzip.open('ref/GCF_036323735.1_GRCr8_protein.faa.gz','rt'),'fasta') if 'schlafen' in r.description.lower()}
rows=[]
for gid,g in G.items():
    best=None
    for mid,par in M.items():
        if par!=gid or mid not in CDS: continue
        pid=CDS[mid]['prot']; L=len(prot.get(pid,''))
        # prefer NP_ over XP_, then longest
        key=(pid.startswith('NP_'),L)
        if best is None or key>best[0]: best=(key,mid,pid)
    _,mid,pid=best; seq=prot[pid].rstrip('*'); pa=ProteinAnalysis(seq.replace('X',''))
    nm=sum(1 for m,p in M.items() if p==gid)
    rows.append(dict(gene=g['gene'],GeneID=g['geneid'],chr=chrname.get(g['seqid'],g['seqid']),start=g['start'],end=g['end'],strand=g['strand'],
        gene_len_bp=g['end']-g['start']+1,n_transcripts=nm,exons=len(EX[mid]),protein=pid,aa=len(seq),MW_kDa=round(pa.molecular_weight()/1000,2),
        pI=round(pa.isoelectric_point(),2),instability=round(pa.instability_index(),1),GRAVY=round(pa.gravy(),3),aliphatic=round(sum(seq.count(x) for x in 'A')/len(seq)*100+2.9*seq.count('V')/len(seq)*100+3.9*(seq.count('I')+seq.count('L'))/len(seq)*100,1)))
    open(f'out/{g["gene"]}.faa','w').write(f'>{g["gene"]}|{pid}\n{seq}\n')
T=pd.DataFrame(rows).sort_values(['chr','start']); T.to_csv('out/slfn_rat_genes.csv',index=False); print(T.to_string(index=False))
