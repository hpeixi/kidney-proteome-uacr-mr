# Circulating plasma proteins and albuminuria (cis-pQTL MR)

Code and derived tables for a two-sample Mendelian randomization study of plasma proteins and albuminuria (UACR), with eGFR as a comparator outcome.

- **Exposure**: SomaScan plasma pQTL, 3,282 proteins, n = 3,301 (OpenGWAS `prot-a`)
- **Outcomes**: CKDGen UACR and eGFR, European ancestry
- **Methods**: cis instruments (±500 kb, P < 5e-8, F > 10), Wald ratio or IVW, approximate-Bayes-factor colocalization, independent replication
- **Results**: 582 proteins with cis instruments; 12 associated with UACR (FDR < 0.05); 8 colocalized (PP.H4 ≥ 0.8); SERPINA10 replicated, BAG3 opposite direction

**Contents**

```
analysis/   8 Python scripts (instruments, MR, colocalization, replication, figures)
results/    11 derived tables
figures/    4 figures (PNG 300 dpi + PDF)
metadata/   CITATION.cff, data dictionary
```

All input summary statistics are public (CKDGen, OpenGWAS, GWAS Catalog) and not redistributed here.

Citation and DOI: see `metadata/CITATION.cff` (DOI pending release).
License: MIT for code, CC BY 4.0 for derived tables and figures.
Contact: hpeixi@foxmail.com
