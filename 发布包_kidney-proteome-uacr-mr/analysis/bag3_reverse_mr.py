import gzip, os
import pandas as pd, numpy as np
from scipy.stats import norm

D = os.path.join(os.path.dirname(__file__), "data")
UACR = os.path.join(D, "UACR_EA_Teumer2019.tbl.rsid.gz")
BAG3 = os.path.join(D, "ukbppp", "GCST90468416.h.tsv.gz")

# 1) UACR 全基因组显著位点
hit = []
with gzip.open(UACR, "rt") as f:
    h = f.readline().split(); j = {c: i for i, c in enumerate(h)}
    for line in f:
        q = line.split()
        if len(q) < len(h):
            continue
        try:
            p = float(q[j["P-value"]])
        except Exception:
            continue
        if p < 5e-8:
            hit.append((q[j["RSID"]].lower(), q[j["Chr"]], int(q[j["Pos_b37"]]),
                        q[j["Allele1"]].upper(), q[j["Allele2"]].upper(),
                        float(q[j["Effect"]]), float(q[j["StdErr"]]), p))
u = pd.DataFrame(hit, columns=["rsid","chr","pos","a1","a2","beta","se","p"])
print("UACR 全基因组显著位点:", len(u))

# 2) 距离剪枝（1 Mb，近似独立）
u = u.sort_values("p")
kept = []
for _, r in u.iterrows():
    if all(not (r["chr"] == k["chr"] and abs(r["pos"] - k["pos"]) < 1_000_000) for k in kept):
        kept.append(r)
iv = pd.DataFrame(kept)
print("1 Mb 剪枝后工具变量:", len(iv))

# 3) 在 BAG3 文件中查这些位点
want = set(iv["rsid"])
rows = []
with gzip.open(BAG3, "rt") as f:
    hh = f.readline().rstrip("\n").split("\t"); k = {c: i for i, c in enumerate(hh)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hh):
            continue
        rs = p[k["rsid"]].lower()
        if rs in want:
            rows.append((rs, p[k["effect_allele"]].upper(), p[k["other_allele"]].upper(),
                         float(p[k["beta"]]), float(p[k["standard_error"]]), float(p[k["p_value"]])))
b = pd.DataFrame(rows, columns=["rsid","ea","nea","beta_b","se_b","p_b"])
print("BAG3 中匹配到:", len(b))

m = iv.merge(b, on="rsid", how="inner")
keep = []
for _, x in m.iterrows():
    if x["a1"] == x["ea"] and x["a2"] == x["nea"]:
        keep.append((x["beta"], x["se"], x["beta_b"], x["se_b"], x["rsid"], x["chr"], x["pos"]))
    elif x["a1"] == x["nea"] and x["a2"] == x["ea"]:
        keep.append((x["beta"], x["se"], -x["beta_b"], x["se_b"], x["rsid"], x["chr"], x["pos"]))
kk = pd.DataFrame(keep, columns=["beta_u","se_u","beta_b","se_b","rsid","chr","pos"])
print("等位基因匹配后:", len(kk))

if len(kk):
    w = (kk.beta_u/kk.se_u)**2
    beta = (w*(kk.beta_b/kk.beta_u)).sum()/w.sum()
    se = np.sqrt((w**2*((kk.se_b/kk.beta_u)**2 + ((kk.beta_b**2)*(kk.se_u**2))/(kk.beta_u**4))).sum())/w.sum()
    p = 2*norm.sf(abs(beta/se))
    out = pd.DataFrame([dict(gene="BAG3", direction="UACR -> BAG3 (reverse)",
                             n_snp=len(kk), beta=beta, se=se, p=p,
                             n_exp=47745)])
    out.to_csv("bag3_reverse_mr.csv", index=False)
    kk.to_csv("bag3_reverse_instruments.csv", index=False)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\n前 5 个工具变量:")
    print(kk.head(5).to_string(index=False, float_format=lambda x: f"{x:.4g}"))
else:
    print("无可用工具变量")
