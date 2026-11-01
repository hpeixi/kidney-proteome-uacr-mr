# Circulating plasma proteins and albuminuria: cis-pQTL Mendelian randomization and colocalization

This repository contains the analysis code, derived result tables and figures for a two-sample
Mendelian randomization (MR) study of circulating plasma proteins and albuminuria.

## Study summary

- Exposure: SomaScan plasma pQTL summary statistics (Sun et al. 2018; 3,301 European-ancestry
  participants; 3,282 aptamer-based protein datasets, accessed through OpenGWAS `prot-a-*`).
- Primary outcome: urinary albumin-to-creatinine ratio (UACR), CKDGen European-ancestry summary
  statistics (Teumer et al. 2019).
- Comparator outcome: eGFR, CKDGen European-ancestry summary statistics (Wuttke et al. 2019).
- Instruments: cis variants within 500 kb of the encoding gene, P < 5e-8, F > 10.
- Estimation: Wald ratio (single instrument) or inverse-variance weighting (>= 2 instruments).
- Colocalization: approximate Bayes factor, +/-100 kb around the lead pQTL, priors 1e-4/1e-4/1e-5.
- Independent replication: SERPINA10 (EBI GCST90019385, n = 10,708), BAG3 (UKB-PPP GCST90468416,
  n = 47,745), ANGPTL1 (EBI GCST90010282, n = 1,313, no genome-wide significant cis variant).

## Key results

| Item | Value |
|---|---|
| Proteins with >= 1 genome-wide significant pQTL | 2,162 / 3,282 |
| Proteins with cis instruments | 582 (769 variant-protein pairs) |
| UACR associations (FDR < 0.05) | 12 |
| eGFR associations (FDR < 0.05) | 31 |
| Colocalized (PP.H4 >= 0.8) | 8 of 12 |
| Replicated (consistent direction) | SERPINA10 |
| Not replicated (opposite direction) | BAG3 |
| Not evaluable | ANGPTL1 |
| No independent data | 5 (APOA5, MCL1, ISLR2, PCOLCE, RSPO3) |

## Repository layout

```
analysis/     analysis scripts (Python 3.12)
results/      derived tables used in the manuscript
figures/      figures (PNG 300 dpi and PDF)
metadata/     data dictionary and citation metadata
```

## Reproducing the analysis

1. Install Python 3.12 with numpy, scipy, pandas, matplotlib, openpyxl.
2. Download the summary statistics listed in `metadata/DATA_DICTIONARY.md` into `data/`.
3. Run the scripts in `analysis/` in the order given in the data dictionary.

The scripts expect the outcome files in `data/` (UACR and eGFR) and the replication files under
`data/ukbppp/`. Paths can be edited at the top of each script.

## Data availability

All summary statistics are publicly available from CKDGen, OpenGWAS, the GWAS Catalog and Synapse;
no individual-level data were used. Derived tables in `results/` are released under CC BY 4.0.

## License

Code: MIT. Derived tables and figures: CC BY 4.0.

## Citation

See `CITATION.cff`. A Zenodo DOI will be minted on the first release; replace the placeholder in
`CITATION.cff` and `metadata/DOI.txt` after publication.
