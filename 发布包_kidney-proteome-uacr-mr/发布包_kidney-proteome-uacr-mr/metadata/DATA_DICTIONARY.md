# Data dictionary

## Inputs (public summary statistics, not redistributed here)

| File | Source | Content |
|---|---|---|
| UACR_EA_Teumer2019.tbl.rsid.gz | CKDGen (Teumer et al. 2019) | UACR GWAS, European ancestry |
| eGFR_EA_Wuttke2019.txt.gz | CKDGen (Wuttke et al. 2019) | eGFR GWAS, European ancestry |
| prot-a-* | OpenGWAS (Sun et al. 2018) | plasma pQTL, 3,282 aptamer datasets |
| GCST90468416.h.tsv.gz | GWAS Catalog (UKB-PPP) | BAG3 pQTL, n = 47,745 |

## Derived tables (results/)

| File | Rows | Columns |
|---|---|---|
| tophits_all.csv | 3,854 | pQTL records per protein-variant |
| outcome_snps.csv | 3,879 | variant effects in UACR and eGFR |
| cis_snps.csv | 769 | cis variant-protein pairs |
| cis_mr_results.csv | 582 | MR estimates per protein |
| candidate_table.csv | 582 | two-outcome classification and direction check |
| 候选蛋白_汇总.csv | 12 | UACR hits with colocalization and replication status |
| coloc_final.csv | 12 | PP.H1-H4 per protein |
| S7_independent_replication.csv | 3 | independent-panel replication results |
| sensitivity_multiinstrument.csv | 145 | IVW restricted to >= 2 instruments |
| bag3_reverse_mr.csv | 1 | reverse MR (UACR -> BAG3) |
| bag3_cis_variants.csv | 405 | BAG3 cis variants from UKB-PPP |

## Scripts (analysis/)

| Script | Purpose |
|---|---|
| cis_mr.py | cis instrument selection and MR |
| coloc_run2.py | approximate Bayes factor colocalization |
| protc_replication.py / independent_replication.py | independent-panel replication |
| bag3_replication.py | BAG3 replication in UKB-PPP |
| bag3_reverse_mr.py | reverse-direction MR for BAG3 |
| sensitivity_multiinstrument.py | multi-instrument sensitivity analysis |
| make_figures.py | figures |
