import json, ssl, urllib.request, urllib.error, os, time, gzip
import pandas as pd, numpy as np
from scipy.stats import norm

TOKEN = open(os.path.join(os.path.dirname(__file__), ".opengwas_token")).read().strip()
ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
BASE = "https://api.opengwas.io/api"
D = os.path.join(os.path.dirname(__file__), "data")

def post(path, payload, timeout=120, retries=3):
    data = json.dumps(payload).encode()
    for a in range(retries):
        try:
            req = urllib.request.Request(BASE+path, data=data, method="POST",
                headers={"Authorization":"Bearer "+TOKEN,"Content-Type":"application/json",
                         "User-Agent":"Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429,500,502,503):
                time.sleep(3*(a+1)); continue
            raise
        except Exception:
            time.sleep(3*(a+1))
    return []

pairs = [("SERPINA10", "ebi-a-GCST90019385"), ("ANGPTL1", "ebi-a-GCST90010282")]
th = []
for _, pid in pairs:
    r = post("/tophits", {"id": [pid]})
    print(pid, "tophits:", len(r))
    th.extend(r)
th = pd.DataFrame(th)
print("rows:", len(th), "cols:", list(th.columns))

want = set(th["rsid"].astype(str).str.lower())
rows = []
with gzip.open(os.path.join(D, "UACR_EA_Teumer2019.tbl.rsid.gz"), "rt") as f:
    header = f.readline().split(); ix = {c: i for i, c in enumerate(header)}
    for line in f:
        p = line.split()
        if len(p) < len(header): continue
        rs = p[ix["RSID"]].lower()
        if rs in want:
            rows.append((rs, p[ix["Chr"]], int(p[ix["Pos_b37"]]), p[ix["Allele1"]].upper(),
                         p[ix["Allele2"]].upper(), float(p[ix["Effect"]]), float(p[ix["StdErr"]])))
out = pd.DataFrame(rows, columns=["rsid","chr","pos_b37","a1","a2","beta","se"])
print("UACR matched:", len(out))

def tss(symbols):
    data = json.dumps({"symbols": symbols}).encode()
    req = urllib.request.Request("https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens",
        data=data, method="POST", headers={"Content-Type":"application/json",
        "Accept":"application/json","User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60, context=ctx) as r:
        return json.loads(r.read().decode())

coords = tss([p[0] for p in pairs])
ref = pd.read_csv("候选蛋白_汇总.csv").set_index("gene")
res = []
for gene, pid in pairs:
    v = coords.get(gene)
    if not isinstance(v, dict):
        res.append(dict(gene=gene, panel=pid, status="无基因坐标")); continue
    chrom = str(v["seq_region_name"]); t = int(v["start"]) if v["strand"] == 1 else int(v["end"])
    h = th[th["id"] == pid].copy()
    h["rsid"] = h["rsid"].astype(str).str.lower()
    h = h.merge(out, on="rsid", how="inner", suffixes=("_p", ""))
    cis = h[(h["chr"] == chrom) & ((h["pos_b37"] - t).abs() <= 500_000)]
    keep = []
    for _, x in cis.iterrows():
        ea, nea = str(x["ea"]).upper(), str(x["nea"]).upper()
        if x["a1"] == ea and x["a2"] == nea: keep.append((float(x["beta_p"]), float(x["se_p"]), x["beta"], x["se"]))
        elif x["a1"] == nea and x["a2"] == ea: keep.append((-float(x["beta_p"]), float(x["se_p"]), x["beta"], x["se"]))
    if not keep:
        res.append(dict(gene=gene, panel=pid, status=f"无可用 cis 位点（候选 {len(cis)}）")); continue
    k = np.array(keep)
    w = (k[:,0]/k[:,1])**2
    b = (w*(k[:,2]/k[:,0])).sum()/w.sum()
    se = np.sqrt((w**2*((k[:,3]/k[:,0])**2 + ((k[:,2]**2)*(k[:,1]**2))/(k[:,0]**4))).sum())/w.sum()
    p = 2*norm.sf(abs(b/se))
    ba = float(ref.loc[gene, "UACR_beta"])
    res.append(dict(gene=gene, panel=pid, n_snp=len(k), beta_rep=b, se_rep=se, p_rep=p,
                    beta_prota=ba, direction="一致" if np.sign(b) == np.sign(ba) else "相反",
                    status="ok"))

r = pd.DataFrame(res)
r.to_csv("independent_replication.csv", index=False)
print(r.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
