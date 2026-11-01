import pandas as pd, gzip, os, json

D = os.path.join(os.path.dirname(__file__), "data")
th = pd.read_csv("tophits_all.csv")
want = set(th["rsid"].astype(str).str.lower())
print("instrument SNPs:", len(want), flush=True)

def scan(fn):
    rows = []
    with gzip.open(os.path.join(D, fn), "rt") as f:
        header = f.readline().split()
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.split()
            if len(p) < len(header):
                continue
            rs = p[idx["RSID"]].lower()
            if rs in want:
                rows.append((rs, p[idx["Chr"]], int(p[idx["Pos_b37"]]),
                             p[idx["Allele1"]].upper(), p[idx["Allele2"]].upper(),
                             float(p[idx["Freq1"]]), float(p[idx["Effect"]]),
                             float(p[idx["StdErr"]]), float(p[idx["P-value"]]),
                             int(float(p[idx["n_total_sum"]]))))
    return pd.DataFrame(rows, columns=["rsid","chr","pos_b37","a1","a2","freq1","beta","se","p","n"])

u = scan("UACR_EA_Teumer2019.tbl.rsid.gz"); u["outcome"] = "UACR"
print("UACR matched SNPs:", len(u), flush=True)
e = scan("eGFR_EA_Wuttke2019.txt.gz");      e["outcome"] = "eGFR"
print("eGFR matched SNPs:", len(e), flush=True)

out = pd.concat([u, e])
out.to_csv("outcome_snps.csv", index=False)
print("saved", len(out), "rows; unique rsid:", out["rsid"].nunique(), flush=True)
print(u.head(3).to_string(index=False), flush=True)
