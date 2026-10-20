import json,requests,time
K=json.load(open('keys.json')); out={}
for k,d in K.items():
    for t in range(4):
        try:
            r=requests.get(f'https://api.crossref.org/works/{d}',timeout=40)
            if r.status_code==200:
                m=r.json()['message']; a=m.get('author',[{}]); 
                out[k]=dict(doi=d,authors=[(x.get('family',''),x.get('given','')) for x in a],title=m.get('title',[''])[0],journal=(m.get('container-title') or [''])[0],year=(m.get('issued',{}).get('date-parts',[[None]])[0][0]),volume=m.get('volume',''),issue=m.get('issue',''),page=m.get('page','') or m.get('article-number',''),type=m.get('type'))
                break
            elif r.status_code==404: out[k]={'doi':d,'error':404}; break
        except Exception as e: time.sleep(2)
    else: out[k]={'doi':d,'error':'fail'}
    time.sleep(0.3)
json.dump(out,open('refs_verified.json','w'),indent=1)
for k,v in out.items(): print(k, v.get('error') or f"{v['authors'][0][0] if v['authors'] else '?'} {v['year']} | {v['title'][:80]} | {v['journal'][:35]}")
