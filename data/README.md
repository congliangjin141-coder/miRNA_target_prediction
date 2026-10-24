# Input data

The input data are not tracked in the repository. Place the files in the folders below.

```
data/
├── raw/            # database downloads and internal expression data, read by 01
├── reference/      # written by 01_generate_input_sets.py, read by 02 and the check
└── predictions/    # the seven tool files, read by 02, the check and Fig2 (when recomputing curves)
```

## `data/raw/`

| File | Content |
|------|---------|
| `Homo_sapiens_TarBase-v9.tsv` | TarBase v9.0 human download, one row per interaction record (4,724,537 records); columns `mirna_name` and `gene_name` are used |
| `hsa_MTI.csv` | miRTarBase 2025 human download; columns `miRNA`, `Target Gene`, `Support Type` are used |
| `bulkRNA_counts_cleaned.csv` | internal iPSC bulk RNA-seq counts, genes × samples (samples `V1_1` … `V10_2`) |
| `miRNAmature_counts.csv` | internal iPSC mature miRNA counts, samples × miRNAs (samples `V1-1` … `V10-2`) |

## `data/reference/` (written by `01`)

| File | Columns | Content |
|------|---------|---------|
| `set1_in_matrix.csv` | `miRNA`, `gene`, `source` | miRTarBase and TarBase pairs inside the expression matrix; `source` is miRTarBase, TarBase or both |
| `set2_A_loose.csv` | `mirna_name`, `gene_name`, `correlation`, `pvalue`, `FDR` | iPSC pairs with Pearson r < −0.2 and BH FDR < 0.05 (256,876 pairs) |
| `set1_inter_set2.csv` | as `set2_A_loose.csv` plus `source` | pairs in both sets (descriptive only) |
| `null_relations_no_interaction.parquet` | `mirna_name`, `gene_name`, `correlation` | pairs with a defined correlation in neither set (57,369,107 pairs) |

## `data/predictions/`

Each file needs the columns `mirna_name`, `gene_name` and `score`; other columns are ignored.
Tools are listed by the number of predicted pairs, fewest to most.

| File | Tool | Stronger prediction |
|------|------|---------------------|
| `TargetScan v8.0.csv` | TargetScan v8.0 | lower score |
| `miRDB v6.0.csv` | miRDB v6.0 | higher score |
| `mirmap.csv` | miRmap | lower score |
| `diana_microT.csv` | DIANA-microT | higher score |
| `rna22.csv` | RNA22 | lower score |
| `mirwalk.csv` | miRWalk | higher score |
| `mirDIP v5.2.csv` | mirDIP v5.2 | higher score |
