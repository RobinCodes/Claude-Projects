import json, urllib.request, urllib.parse, time, os, hashlib

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'oeis_cache')
os.makedirs(CACHE, exist_ok=True)

def query(terms, maxn=10):
    """terms: list of ints. Returns list of (Anum, name)."""
    t = [str(x) for x in terms[:maxn]]
    q = ",".join(t)
    h = hashlib.md5(q.encode()).hexdigest()
    fp = os.path.join(CACHE, h+'.json')
    if os.path.exists(fp):
        return json.load(open(fp))
    url = "https://oeis.org/search?q=" + urllib.parse.quote(q) + "&fmt=json"
    req = urllib.request.Request(url, headers={'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) research/1.0'})
    try:
        with urllib.request.urlopen(req, timeout=25) as r:
            data = json.load(r)
    except Exception as e:
        return [('ERR', str(e))]
    res = []
    if isinstance(data, list):
        for d in data[:6]:
            res.append(('A%06d' % d['number'], d.get('name','')))
    json.dump(res, open(fp,'w'))
    time.sleep(0.6)
    return res

if __name__ == '__main__':
    import sys
    print(query([int(x) for x in sys.argv[1].split(',')]))
