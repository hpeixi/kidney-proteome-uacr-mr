import json, ssl, urllib.request, urllib.parse, urllib.error, os, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
CACHE = "gene_map.json"
lock = threading.Lock()
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}

th = pd.read_csv("tophits_all.csv")
info = json.load(open("ogwas_all.json"))
ids = sorted(th["id"].unique())
print("proteins with hits:", len(ids), "cached:", len(cache), flush=True)

def query(name):
    for reviewed in (True, False):
        sp = " AND organism_id:9606" + (" AND reviewed:true" if reviewed else "")
        q = urllib.parse.quote(f'protein_name:"{name}"{sp}')
        url = ("https://rest.uniprot.org/uniprotkb/search?query=" + q +
               "&fields=gene_primary&format=tsv&size=1")
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=12, context=ctx) as r:
                txt = r.read().decode("utf-8","replace").strip().splitlines()
            if len(txt) >= 2:
                g = txt[1].split("\t")[0]
                if g and g != "-":
                    return g
        except Exception:
            pass
    return None

todo = [k for k in ids if k not in cache or not cache[k].get("gene")]
print("todo:", len(todo), flush=True)

def work(k):
    return k, query(info[k]["trait"])

t0 = time.time(); done = 0
with ThreadPoolExecutor(max_workers=16) as ex:
    futs = {ex.submit(work, k): k for k in todo}
    for f in as_completed(futs):
        k, g = f.result()
        with lock:
            cache[k] = {"gene": g}
            done += 1
            if done % 200 == 0:
                json.dump(cache, open(CACHE, "w"))
                el = time.time()-t0
                print(f"{done}/{len(todo)} {done/el:.1f}/s ETA {(len(todo)-done)/(done/el)/60:.1f} min", flush=True)

json.dump(cache, open(CACHE, "w"))
ok = sum(1 for k in ids if cache.get(k, {}).get("gene"))
print("FINISHED mapped", ok, "/", len(ids), flush=True)
