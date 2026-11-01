import gzip, os, ssl, urllib.request
import pandas as pd, numpy as np
from scipy.stats import norm

D = os.path.join(os.path.dirname(__file__), "data")
F = os.path.join(D, "ukbppp", "GCST90468416.h.tsv.gz")
UACR = os.path.join(D, "UACR_EA_Teumer2019.tbl.rsid.gz")
BAG3_TSS_B37 = 121_410_882   # Ensembl GRCh37, chr10
WIN = 500_000

print("BAG3 file:", os.path.getsize(F) / 1e6, "MB")

# 1) 从 UKB-PPP 提取 chr10 全基因组显著位点
rows = []
with gzip.open(F, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(header)}
    print("columns:", header)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(header):
            continue
        if p[ix["chromosome"]] != "10":
            continue
        try:
            if float(p[ix["p_value"]]) >= 5e-8:
                continue
        except Exception:
            continue
        rs = p[ix["rsid"]]
        if rs in ("NA", "", "."):
            continue
        rows.append((rs, float(p[ix["beta"]]), float(p[ix["standard_error"]]),
                     p[ix["effect_allele"]].upper(), p[ix["other_allele"]].upper(),
                     float(p[ix["effect_allele_frequency"]]) if p[ix["effect_allele_frequency"]] not in ("NA", "") else np.nan,
                     float(p[ix["p_value"]]), int(float(p[ix["n"]])) if p[ix["n"]] not in ("NA", "") else np.nan))
ukb = pd.DataFrame(rows, columns=["rsid","beta_p","se_p","ea","nea","eaf","p","n"])
print("chr10 genome-wide significant with rsid:", len(ukb))
ukb.to_csv("bag3_ukbppp_chr10_sig.csv", index=False)

# 2) 在 UACR 中查这些位点（GRCh37 坐标）
want = set(ukb["rsid"].str.lower())
u = []
with gzip.open(UACR, "rt") as f:
    h = f.readline().split(); j = {c: i for i, c in enumerate(h)}
    for line in f:
        q = line.split()
        if len(q) < len(h):
            continue
        rs = q[j["RSID"]].lower()
        if rs in want:
            u.append((rs, q[j["Chr"]], int(q[j["Pos_b37"]]), q[j["Allele1"]].upper(),
                      q[j["Allele2"]].upper(), float(q[j["Effect"]]), float(q[j["StdErr"]])))
uacr = pd.DataFrame(u, columns=["rsid","chr","pos_b37","a1","a2","beta_u","se_u"])
print("matched in UACR:", len(uacr))

# 3) cis 限定（b37）
cis = uacr[(uacr["chr"] == "10") & ((uacr["pos_b37"] - BAG3_TSS_B37).abs() <= WIN)]
m = cis.merge(ukb.assign(rsid=ukb.rsid.str.lower()), on="rsid", how="inner")
print("cis variants:", len(m))

keep = []
for _, x in m.iterrows():
    if x["a1"] == x["ea"] and x["a2"] == x["nea"]:
        keep.append((x["beta_p"], x["se_p"], x["beta_u"], x["se_u"], x["rsid"], x["pos_b37"]))
    elif x["a1"] == x["nea"] and x["a2"] == x["ea"]:
        keep.append((-x["beta_p"], x["se_p"], x["beta_u"], x["se_u"], x["rsid"], x["pos_b37"]))
k = pd.DataFrame(keep, columns=["beta_p","se_p","beta_u","se_u","rsid","pos_b37"])
print("harmonised:", len(k))
if len(k):
    k["F"] = (k.beta_p/k.se_p)**2
    w = (k.beta_p/k.se_p)**2
    b = (w*(k.beta_u/k.beta_p)).sum()/w.sum()
    se = np.sqrt((w**2*((k.se_u/k.beta_p)**2 + ((k.beta_u**2)*(k.se_p**2))/(k.beta_p**4))).sum())/w.sum()
    p = 2*norm.sf(abs(b/se))
    res = pd.DataFrame([dict(gene="BAG3", panel="UKB-PPP GCST90468416", n_snp=len(k),
                             lead=k.sort_values("F", ascending=False).rsid.iloc[0],
                             beta_rep=b, se_rep=se, p_rep=p,
                             beta_prota=0.05145, se_prota=0.01545,
                             direction="一致" if np.sign(b) == np.sign(0.05145) else "相反")])
    res.to_csv("bag3_replication.csv", index=False)
    k.to_csv("bag3_cis_variants.csv", index=False)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\ncis variants detail:")
    print(k.sort_values("F", ascending=False).head(10).to_string(index=False, float_format=lambda x: f"{x:.4g}"))
else:
    print("no harmonised cis variant")
