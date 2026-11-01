import json, ssl, urllib.request, urllib.error, os, sys, time
import pandas as pd, numpy as np
from scipy.stats import norm, chi2

ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
GRCH37 = "https://grch37.rest.ensembl.org"

def post_symbols(symbols, timeout=90):
    data = json.dumps({"symbols": symbols}).encode()
    req = urllib.request.Request(GRCH37 + "/lookup/symbol/homo_sapiens", data=data,
        method="POST", headers={"Content-Type":"application/json",
                                "Accept":"application/json","User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
        return json.loads(r.read().decode())

def gene_tss(symbols):
    out = {}
    for i in range(0, len(symbols), 180):
        chunk = symbols[i:i+180]
        try:
            res = post_symbols(chunk)
        except Exception as e:
            print("  ensembl chunk fail", i, repr(e)[:120], flush=True); continue
        if isinstance(res, dict):
            for s, v in res.items():
                if isinstance(v, dict):
                    out[s] = (str(v["seq_region_name"]), int(v["start"]) if v["strand"] == 1 else int(v["end"]))
        time.sleep(0.3)
    return out

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "test":
    print(gene_tss(["APOE", "UMOD", "FGF5", "FCGR2A"]))
    sys.exit()

MAXKB = 500_000
gm = json.load(open("gene_map.json"))
th = pd.read_csv("tophits_all.csv")

sym2id = {}
for pid, v in gm.items():
    g = (v or {}).get("gene")
    if g:
        sym2id.setdefault(g, []).append(pid)

print("genes to locate:", len(sym2id), flush=True)
tss = gene_tss(sorted(sym2id))
print("located:", len(tss), flush=True)

oc = pd.read_csv("outcome_snps.csv")
oc = oc.drop_duplicates(subset=["rsid", "outcome"])
piv = {o: g.set_index("rsid") for o, g in oc.groupby("outcome")}
b37 = oc.drop_duplicates("rsid")[["rsid", "chr", "pos_b37"]].copy()
b37["chr"] = b37["chr"].astype(str)
b37["rsid"] = b37["rsid"].str.lower()

rows = []
for gene, pids in sym2id.items():
    if gene not in tss:
        continue
    chrom, t = tss[gene]
    for pid in pids:
        hits = th[th["id"] == pid].copy()
        hits["rsid"] = hits["rsid"].astype(str).str.lower()
        hits = hits.merge(b37, on="rsid", how="inner", suffixes=("", "_b37"))
        cis = hits[(hits["chr_b37"] == chrom) & ((hits["pos_b37"] - t).abs() <= MAXKB)]
        for _, h in cis.iterrows():
            rs = h["rsid"]
            rec = dict(gene=gene, id=pid, protein=h["trait"], rsid=rs,
                       b_exp=float(h["beta"]), se_exp=float(h["se"]),
                       ea=h["ea"].upper(), nea=h["nea"].upper())
            for o in ("UACR", "eGFR"):
                tab = piv.get(o)
                if tab is None or rs not in tab.index:
                    rec[o+"_b"] = np.nan; rec[o+"_se"] = np.nan; continue
                r = tab.loc[rs]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                bo, so = float(r["beta"]), float(r["se"])
                if r["a1"] == rec["ea"] and r["a2"] == rec["nea"]:
                    pass
                elif r["a1"] == rec["nea"] and r["a2"] == rec["ea"]:
                    bo = -bo
                else:
                    rec[o+"_b"] = np.nan; rec[o+"_se"] = np.nan; continue
                rec[o+"_b"] = bo; rec[o+"_se"] = so
            rows.append(rec)

snps = pd.DataFrame(rows)
print("cis SNP-protein pairs:", len(snps), flush=True)
snps.to_csv("cis_snps.csv", index=False)

res = []
for (pid, gene), g in snps.groupby(["id", "gene"]):
    rec = dict(id=pid, gene=gene, protein=g["protein"].iloc[0], n_snp=len(g))
    F = (g["b_exp"]/g["se_exp"])**2
    rec["F_mean"] = F.mean(); rec["F_min"] = F.min()
    for o in ("UACR", "eGFR"):
        sub = g.dropna(subset=[o+"_b", o+"_se"])
        if len(sub) == 0:
            rec[o+"_beta"] = np.nan; rec[o+"_se"] = np.nan; rec[o+"_p"] = np.nan; continue
        w = (sub["b_exp"]/sub["se_exp"])**2
        b = (w * (sub[o+"_b"]/sub["b_exp"])).sum()/w.sum()
        se = np.sqrt((w**2 * ((sub[o+"_se"]/sub["b_exp"])**2 + ((sub[o+"_b"]**2)*(sub["se_exp"]**2))/(sub["b_exp"]**4))).sum())/w.sum()
        rec[o+"_beta"] = b; rec[o+"_se"] = se; rec[o+"_p"] = 2*norm.sf(abs(b/se))
    res.append(rec)

rd = pd.DataFrame(res)
def bh(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p); q = np.empty(n); prev = 1
    for rank in range(n-1, -1, -1):
        i = order[rank]; prev = min(prev, p[i]*n/(rank+1)); q[i] = prev
    return q

for o in ("UACR", "eGFR"):
    m = rd[o+"_p"].notna()
    rd.loc[m, o+"_fdr"] = bh(rd.loc[m, o+"_p"])

rd.to_csv("cis_mr_results.csv", index=False)
print("proteins with cis instruments:", len(rd), flush=True)
print("UACR FDR<0.05:", int((rd["UACR_fdr"] < 0.05).sum()),
      "eGFR FDR<0.05:", int((rd["eGFR_fdr"] < 0.05).sum()), flush=True)
cols = ["protein","gene","n_snp","F_min","UACR_beta","UACR_se","UACR_p","UACR_fdr"]
print(rd.nsmallest(15, "UACR_p")[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
