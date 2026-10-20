import requests,json,subprocess,pandas as pd
from Bio import SeqIO
s=requests.Session(); H={'Content-Type':'application/json'}
ref=json.load(open('out/ref_pfam.json'))
seqs={}
for k,v in ref.items():
    sp,sym=k.split('_',2)[0]+'_'+k.split('_',2)[1],k.split('_',2)[2]
    tag={'homo_sapiens':'Hsa','mus_musculus':'Mmu'}[sp]; seqs[f'{tag}_{sym}']=v['seq']
for g in pd.read_csv('out/slfn_rat_genes.csv').gene:
    seqs[f'Rno_{g}']=str(next(SeqIO.parse(f'out/{g}.faa','fasta')).seq)
# extra species: dog, cattle via human orthologues (canonical)
for sym in ['SLFN5','SLFN11','SLFN12','SLFN12L','SLFN13','SLFN14','SLFNL1']:
    for sp,tag in [('canis_lupus_familiaris','Cfa'),('bos_taurus','Bta')]:
        h=s.get(f'https://rest.ensembl.org/homology/symbol/homo_sapiens/{sym}',headers=H,params={'target_species':sp,'type':'orthologues','sequence':'protein','aligned':0}).json()
        try:
            for d in h['data'][0]['homologies']:
                t=d['target']; nm=f"{tag}_{sym}like_{t['protein_id'][-6:]}"
                seqs[nm]=t['seq']
        except Exception as e: print('err',sym,sp,str(h)[:120])
open('out/slfn_tree_input.fa','w').write(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()))
print(len(seqs)); print(sorted(seqs))
