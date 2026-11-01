import json, ssl, urllib.request, urllib.error, os, time
import pandas as pd, numpy as np
from scipy.special import logsumexp

TOKEN = open(os.path.join(os.path.dirname(__file__), ".opengwas_token")).read().strip()
ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
BASE = "https://api.opengwas.io/api"
WIN = 100_000

def post(payload, timeout=120, retries=3):
    data = json.dumps(payload).encode()
    last = None
    for a in range(retries):
        try:
            req = urllib.request.Request(BASE+"/associations", data=data, method="POST",
                headers={"Authorization":"Bearer "+TOKEN,"Content-Type":"application/json",
                         "User-Agent":"Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (429,500,502,503):
                time.sleep(2*(a+1)); continue
            raise
        except Exception as e:
            last = e; time.sleep(2*(a+1))
    if last: raise last
    return []

def fetch(pid, variants):
    out, stack = [], [variants]
    while stack:
        ch = stack.pop()
        if not ch: continue
        try:
            out.extend(post({"id": [pid], "variant": ch}, retries=2))
        except urllib.error.HTTPError as e:
            if e.code == 400 and len(ch) > 1:
                m = len(ch)//2; stack.extend([ch[:m], ch[m:]])
            elif e.code == 400:
                continue
        time.sleep(0.15)
    return out

def abf(b, s):
    V = s**2; z = b/s; W = np.maximum(b**2 - s**2, 1e-8)
    return 0.5*(np.log(V) - np.log(V+W) + (z**2)*W/(V+W))

def coloc(b1, s1, b2, s2):
    l1, l2 = abf(b1, s1), abf(b2, s2)
    a, b, ab = logsumexp(l1), logsumexp(l2), logsumexp(l1+l2)
    p1 = p2 = 1e-4; p12 = 1e-5
    lh3 = (a+b) + np.log1p(-np.exp(ab-(a+b))) if a+b > ab else -np.inf
    lp = np.array([np.log(p1)+np.log(p2), np.log(p1)+a, np.log(p2)+b,
                   np.log(p1)+np.log(p2)+lh3, np.log(p12)+ab])
    p = np.exp(lp - logsumexp(lp))
    return p

cis = pd.read_csv("cis_mr_results.csv")
cis = cis[(cis["UACR_fdr"] < 0.05)].copy()
reg = pd.read_csv("uacr_regions.csv")
snps = pd.read_csv("cis_snps.csv")
pos = pd.read_csv("outcome_snps.csv").drop_duplicates("rsid")[["rsid","chr","pos_b37"]]

out = []
for _, r in cis.iterrows():
    pid, gene = r["id"], r["gene"]
    s = snps[snps["id"] == pid].copy()
    if s.empty:
        continue
    s["z2"] = (s["b_exp"]/s["se_exp"])**2
    lead = s.sort_values("z2", ascending=False).iloc[0]
    lead_rs = lead["rsid"]
    p = pos[pos["rsid"] == lead_rs]
    if p.empty:
        print(gene, "lead not in outcome file"); continue
    chrom, lpos = str(p["chr"].iloc[0]), int(p["pos_b37"].iloc[0])
    sub = reg[(reg["gene"] == gene) & (reg["chr"].astype(str) == chrom) & ((reg["pos_b37"] - lpos).abs() <= WIN)]
    if len(sub) < 10:
        print(gene, "too few region SNPs", len(sub)); continue
    variants = [f"{chrom}:{int(x)}" for x in sub["pos_b37"]]
    rows = []
    for i in range(0, len(variants), 100):
        rows.extend(fetch(pid, variants[i:i+100]))
    if not rows:
        print(gene, "no pQTL rows"); continue
    q = pd.DataFrame(rows)
    q["rsid"] = q["rsid"].str.lower()
    q = q.rename(columns={"beta":"b1","se":"s1","ea":"ea1","nea":"nea1"})
    m = sub.merge(q[["rsid","b1","s1","ea1","nea1"]], on="rsid", how="inner")
    keep = []
    for _, x in m.iterrows():
        if x["a1"] == x["ea1"].upper() and x["a2"] == x["nea1"].upper():
            keep.append((x["b1"], x["s1"], x["beta"], x["se"]))
        elif x["a1"] == x["nea1"].upper() and x["a2"] == x["ea1"].upper():
            keep.append((-x["b1"], x["s1"], x["beta"], x["se"]))
    if len(keep) < 10:
        print(gene, "too few harmonised", len(keep)); continue
    k = np.array(keep)
    pp = coloc(k[:,0], k[:,1], k[:,2], k[:,3])
    rec = dict(gene=gene, id=pid, n_snp=len(k), lead=lead_rs,
               H1=pp[1], H2=pp[2], H3=pp[3], H4=pp[4])
    out.append(rec)
    print(f"{gene}: n={len(k)} H3={pp[3]:.3f} H4={pp[4]:.3f}", flush=True)

pd.DataFrame(out).to_csv("coloc_final.csv", index=False)
print("saved coloc_final.csv", len(out))
