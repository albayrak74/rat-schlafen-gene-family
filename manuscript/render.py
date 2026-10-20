import json,re,html,unicodedata
R=json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'refs_verified.json')))
def clean(s):
    s=html.unescape(re.sub(r'<[^>]+>','',s)); s=re.sub(r'\s+',' ',s).strip(); return s.replace('tRNA/rRNA- targeting','tRNA/rRNA-targeting')
def initials(g):
    out=''
    for p in re.split(r'\s+',g.strip()):
        if not p: continue
        out+='-'.join(x[0]+'.' for x in p.split('-') if x)
    return out
def fixname(f): return f.replace('‐','-')
info={}
for k,r in R.items():
    au=[a for a in r['authors'] if not a[0].lower().startswith('on behalf')]
    m=re.search(r'(\d{4})',k); y=int(m.group(1)) if m else r['year']
    if k=='Bokobza2025': y=2026
    fam=[fixname(a[0]) for a in au]
    cite=fam[0] if len(fam)==1 else (f'{fam[0]} and {fam[1]}' if len(fam)==2 else f'{fam[0]} et al.')
    names=[f'{fixname(a[0])}, {initials(a[1])}' if a[1] else fixname(a[0]) for a in au]
    auth=', '.join(names[:6])+', et al.' if len(names)>10 else ', '.join(names)
    if k=='Bokobza2025': auth+=', on behalf of the PREMSTEM consortium'
    t=clean(r['title']); j=clean(r['journal']); pg=str(r['page']).replace('-','–'); vol=r.get('volume','')
    if r.get('type')=='book-chapter': ref=f"{auth}, {y}. {t}. {j}, pp. {pg}. https://doi.org/{r['doi']}"
    else: ref=f"{auth}, {y}. {t}{'' if t.endswith(('.','?')) else '.'} {j} {vol}, {pg}. https://doi.org/{r['doi']}"
    sk=''.join(c for c in unicodedata.normalize('NFKD',fam[0]).lower() if not unicodedata.combining(c)).replace(' ','')
    info[k]=dict(cite=cite,year=y,ref=ref,sortkey=(sk,min(len(fam),3),fam[1].lower() if len(fam)==2 else '',y,t))
# disambiguate same cite+year
from collections import defaultdict
g=defaultdict(list)
for k,v in info.items(): g[(v['cite'],v['year'])].append(k)
for ks in g.values():
    if len(ks)>1:
        for i,k in enumerate(sorted(ks,key=lambda k:info[k]['sortkey'])):
            info[k]['ysuf']=str(info[k]['year'])+'abcdef'[i]
for v in info.values():
    v.setdefault('ysuf',str(v['year'])); v['ref']=v['ref'].replace(f", {v['year']}. ",f", {v['ysuf']}. ",1)
def paren(keys):
    ks=sorted(keys,key=lambda k:(info[k]['year'],info[k]['sortkey']))
    return '('+'; '.join(f"{info[k]['cite']}, {info[k]['ysuf']}" for k in ks)+')'
def render(t):
    t=re.sub(r'(\w[\w\s.]*?) \[-@([A-Za-z0-9]+)\]',lambda m: f"{m.group(1)} ({info[m.group(2)]['ysuf']})",t)
    t=re.sub(r'\[((?:@[A-Za-z0-9]+;?\s*)+)\]',lambda m: paren(re.findall(r'@([A-Za-z0-9]+)',m.group(1))),t)
    return t
def reflist(): return [info[k]['ref'] for k in sorted(info,key=lambda k:info[k]['sortkey'])]
if __name__=='__main__':
    for r in reflist(): print(r)
