import pandas as pd, numpy as np
from scipy.stats import norm

s = pd.read_csv("cis_snps.csv")
print("cis_snps:", s.shape, "| 列:", list(s.columns))

def ivw(g):
    w = (g["b_exp"]/g["se_exp"])**2
    b = (w*(g["UACR_b"]/g["b_exp"])).sum()/w.sum()
    se = np.sqrt((w**2*((g["UACR_se"]/g["b_exp"])**2 +
                       ((g["UACR_b"]**2)*(g["se_exp"]**2))/(g["b_exp"]**4))).sum())/w.sum()
    return b, se

rows = []
for (pid, gene), g in s.groupby(["id", "gene"]):
    g = g.dropna(subset=["UACR_b", "UACR_se"])
    if len(g) < 2:
        continue
    b, se = ivw(g)
    rows.append(dict(gene=gene, id=pid, protein=g["protein"].iloc[0], n_snp=len(g),
                     beta=b, se=se, p=2*norm.sf(abs(b/se))))

d = pd.DataFrame(rows)
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for r in range(n-1, -1, -1):
        i = o[r]; prev = min(prev, p[i]*n/(r+1)); q[i] = prev
    return q
d["fdr"] = bh(d["p"])
d = d.sort_values("p")
d.to_csv("sensitivity_multiinstrument.csv", index=False)

full = pd.read_csv("cis_mr_results.csv")
full_hits = set(full[full.UACR_fdr < 0.05].gene)
print("多工具变量蛋白:", len(d), "| FDR<0.05:", int((d.fdr < 0.05).sum()))
print("12 个命中中属于多工具变量的:", sorted(full_hits & set(d.gene)))
print("\n多工具变量子集里 FDR<0.05 的蛋白:")
print(d[d.fdr < 0.05][["gene","n_snp","beta","se","p","fdr"]].to_string(index=False,
      float_format=lambda x: f"{x:.4g}"))
print("\n原 12 个命中在该子集中的结果:")
print(d[d.gene.isin(full_hits)][["gene","n_snp","beta","p","fdr"]].to_string(index=False,
      float_format=lambda x: f"{x:.4g}"))
