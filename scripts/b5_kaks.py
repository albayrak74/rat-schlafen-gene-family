import gzip,re,json,subprocess,pandas as pd,itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import codonalign
G=pd.read_csv('out/slfn_rat_genes.csv'); want=set(G.protein)
rc={}
for r in SeqIO.parse(gzip.open('ref/GCF_036323735.1_GRCr8_cds_from_genomic.fna.gz','rt'),'fasta'):
    m=re.search(r'protein_id=([^\]]+)',r.description)
    if m and m.group(1) in want: rc[m.group(1)]=str(r.seq)
cds={g.gene:rc[g.protein] for _,g in G.iterrows()}
ref=json.load(open('out/ref_pfam.json')); mpid={k.split('_',2)[2]:v['pid'] for k,v in ref.items() if k.startswith('mus')}
mc={}
for r in SeqIO.parse(gzip.open('ref/Mus_musculus.GRCm39.cds.all.fa.gz','rt'),'fasta'):
    m=re.search(r'transcript:(\S+)',r.description); 
    for sym,pid in mpid.items():
        pass
# map mouse translation id -> cds via pep header transcript
pep2tr={}
for r in SeqIO.parse(gzip.open('ref/Mus_musculus.GRCm39.pep.all.fa.gz','rt'),'fasta'):
    pid=r.id.split('.')[0]
    if pid in mpid.values(): pep2tr[pid]=re.search(r'transcript:(\S+)',r.description).group(1).split('.')[0]
for r in SeqIO.parse(gzip.open('ref/Mus_musculus.GRCm39.cds.all.fa.gz','rt'),'fasta'):
    t=r.id.split('.')[0]
    for sym,pid in mpid.items():
        if pep2tr.get(pid)==t: mc['Mmu_'+sym]=str(r.seq)
def kaks(n1,s1,n2,s2):
    s1=s1[:len(s1)//3*3]; s2=s2[:len(s2)//3*3]
    p1=str(Seq(s1).translate()).rstrip('*'); p2=str(Seq(s2).translate()).rstrip('*')
    s1=s1[:len(p1)*3]; s2=s2[:len(p2)*3]
    open('rev/p.fa','w').write(f'>a\n{p1}\n>b\n{p2}\n')
    aln=subprocess.run('mafft --quiet --localpair --maxiterate 1000 rev/p.fa',shell=True,capture_output=True,text=True).stdout
    recs=list(SeqIO.parse(__import__('io').StringIO(aln),'fasta'))
    pa=MultipleSeqAlignment([SeqRecord(r.seq,id=r.id) for r in recs])
    ca=codonalign.build(pa,[SeqRecord(Seq(s1),id='a'),SeqRecord(Seq(s2),id='b')])
    dn,ds=codonalign.codonseq.cal_dn_ds(ca[0],ca[1],method='NG86')
    return round(dn,4),round(ds,4),round(dn/ds,3) if ds>0 else float('nan')
pairs=[('Slfn1','Mmu_Slfn1'),('Slfn2','Mmu_Slfn2'),('Slfn3','Mmu_Slfn3'),('Slfn4','Mmu_Slfn4'),('Slfn5','Mmu_Slfn5'),('Slfn9','Mmu_Slfn9'),('Slfn9','Mmu_Slfn8'),('Slfn14','Mmu_Slfn14'),('Slfnl1','Mmu_Slfnl1')]
rows=[]
for a,b in pairs:
    if b in mc: rows.append(['orthologue','Rno_'+a,b,*kaks(a,cds[a],b,mc[b])])
for a,b in [('Slfn3','Slfn4'),('Slfn1','Slfn2'),('Slfn5','Slfn9'),('Slfn9','Slfn14'),('Slfn5','Slfn14')]:
    rows.append(['paralogue','Rno_'+a,'Rno_'+b,*kaks(a,cds[a],b,cds[b])])
for a,b in [('Mmu_Slfn3','Mmu_Slfn4'),('Mmu_Slfn8','Mmu_Slfn9')]:
    rows.append(['mouse paralogue',a,b,*kaks(a,mc[a],b,mc[b])])
for a,b in [('Slfn3','Mmu_Slfn4'),('Slfn4','Mmu_Slfn3')]:
    rows.append(['cross pair','Rno_'+a,b,*kaks(a,cds[a],b,mc[b])])
T=pd.DataFrame(rows,columns=['type','gene1','gene2','Ka','Ks','Ka_Ks']); T.to_csv('out/slfn_kaks.csv',index=False); print(T.to_string(index=False))
