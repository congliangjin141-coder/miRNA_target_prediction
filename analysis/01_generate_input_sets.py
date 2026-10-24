# ==============================================================================
# 01_generate_input_sets.py - reference sets and candidate pairs
# ==============================================================================
# Author : Congliang Jin
# Date   : 2026-10-23
#
# Description:
#   The three reference sources are read from their raw files. miRTarBase and TarBase
#   list one row per interaction record, so both are collapsed to unique miRNA-gene pairs.
#   set1  miRTarBase (Non-Functional MTI and non-human records removed) and TarBase
#         (all records) pairs inside the iPSC expression matrix;
#         source = miRTarBase, TarBase or both
#   set2  iPSC pairs: Pearson correlation of raw counts over 20 samples (10 donors x
#         2 replicates), r < -0.2 and BH FDR < 0.05 (BH applied to the r < -0.2 pairs)
#   null  pairs with a defined correlation that are in neither set (reference score 0)
#
# Input:
#   data/raw/hsa_MTI.csv                     miRTarBase 2025
#   data/raw/Homo_sapiens_TarBase-v9.tsv     TarBase v9.0
#   data/raw/bulkRNA_counts_cleaned.csv      iPSC RNA-seq counts (genes x samples)
#   data/raw/miRNAmature_counts.csv          iPSC miRNA counts (samples x miRNAs)
#
# Output (data/reference/):
#   set1_in_matrix.csv
#   set2_A_loose.csv
#   set1_inter_set2.csv  (descriptive only)
#   null_relations_no_interaction.parquet
#
# Memory:
#   large: all genes x all miRNAs correlation
#
# Run from the repository root:
#   python analysis/01_generate_input_sets.py
# ==============================================================================
import csv
import gc
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.stats import false_discovery_control
from scipy.stats import t as tdist

PATH_MTB = 'data/raw/hsa_MTI.csv'                  # miRTarBase download
PATH_TB = 'data/raw/Homo_sapiens_TarBase-v9.tsv'   # TarBase download
PATH_RNA = 'data/raw/bulkRNA_counts_cleaned.csv'   # genes x samples
PATH_MIR = 'data/raw/miRNAmature_counts.csv'       # samples x miRNAs
OUT = 'data/reference/'

R_CUT, FDR_CUT = -0.2, 0.05
SAMPLES = [f'V{i}_{j}' for i in range(1, 11) for j in (1, 2)]   # 10 donors x 2 replicates
CHUNK = 3_000_000
os.makedirs(OUT, exist_ok=True)


def clean_mi(s):
    """Strip spaces and lower-case the species prefix (Hsa- -> hsa-)."""
    s = s.strip()
    if '-' in s:
        pre, rest = s.split('-', 1)
        return pre.lower() + '-' + rest
    return s


print('Loading expression matrices ...')
mir_df = pd.read_csv(PATH_MIR, index_col=0)
mir_df.index = mir_df.index.str.replace('-', '_')   # V1-1 -> V1_1
bulk = pd.read_csv(PATH_RNA, index_col=0)
M, G = set(mir_df.columns), set(bulk.index)
print(f'  matrix: {len(M):,} miRNAs x {len(G):,} genes')

# ---------------------------------------------------------------- set1: miRTarBase and TarBase
print('\n[1/4] set1_in_matrix.csv')
mtb = pd.read_csv(PATH_MTB, usecols=['miRNA', 'Target Gene', 'Support Type'], dtype=str)
mtb = mtb[~mtb['Support Type'].str.startswith('Non-Functional')]
mtb['miRNA'] = mtb['miRNA'].apply(clean_mi)
mtb['Target Gene'] = mtb['Target Gene'].str.strip()
mtb = mtb[mtb['miRNA'].str.startswith('hsa-')]
mtb_pairs = set(zip(mtb['miRNA'], mtb['Target Gene']))

tb_pairs, n_tb = set(), 0
with open(PATH_TB, 'r', encoding='utf-8-sig', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    header = next(r)
    mi_col, ge_col = header.index('mirna_name'), header.index('gene_name')
    for row in r:
        n_tb += 1
        if len(row) <= max(mi_col, ge_col):
            continue
        mi, ge = row[mi_col].strip(), row[ge_col].strip()
        if mi and ge:
            tb_pairs.add((clean_mi(mi), ge))
print(f'  miRTarBase pairs: {len(mtb_pairs):,}   TarBase records: {n_tb:,}, pairs: {len(tb_pairs):,}')

rows = [(mi, ge, 'both' if (mi, ge) in tb_pairs else 'miRTarBase')
        for mi, ge in mtb_pairs if mi in M and ge in G]
rows += [(mi, ge, 'TarBase') for mi, ge in tb_pairs
         if (mi, ge) not in mtb_pairs and mi in M and ge in G]
s1m = pd.DataFrame(rows, columns=['miRNA', 'gene', 'source'])
s1m.to_csv(OUT + 'set1_in_matrix.csv', index=False)
print(f'  {len(s1m):,} pairs')

# ---------------------------------------------------------------- set2: iPSC negative correlation
# Pearson correlation of raw counts over the 20 samples, all genes x all miRNAs;
# BH correction applied to the pairs with r < -0.2.
print('\n[2/4] set2_A_loose.csv')
rna, mir = bulk.T.loc[SAMPLES], mir_df.loc[SAMPLES]
n = len(SAMPLES)
X = rna.values.astype(np.float32); Y = mir.values.astype(np.float32)
Xc = X - X.mean(0); Yc = Y - Y.mean(0)
with np.errstate(divide='ignore', invalid='ignore'):
    R = (Xc.T @ Yc) / np.outer(np.sqrt((Xc ** 2).sum(0)), np.sqrt((Yc ** 2).sum(0)))
del X, Y, Xc, Yc; gc.collect()

G_arr, M_arr = rna.columns.to_numpy(), mir.columns.to_numpy()
nM = len(M_arr)
rflat = R.ravel(); del R; gc.collect()
finite = np.isfinite(rflat)            # pairs with a defined correlation = evaluable pairs
rr = np.clip(rflat, -1 + 1e-7, 1 - 1e-7)
with np.errstate(divide='ignore', invalid='ignore'):
    tstat = rr * np.sqrt((n - 2) / (1 - rr ** 2))
p = np.full_like(rflat, np.nan)
p[finite] = 2 * tdist.sf(np.abs(tstat[finite]), n - 2)
del rr, tstat; gc.collect()

neg = finite & (rflat < R_CUT)
pos = np.flatnonzero(neg)
p_neg = p[neg]
s2 = pd.DataFrame({'mirna_name': M_arr[pos % nM], 'gene_name': G_arr[pos // nM],
                   'correlation': rflat[neg], 'pvalue': p_neg,
                   'FDR': false_discovery_control(p_neg, method='bh')})
s2 = s2[s2['FDR'] < FDR_CUT].sort_values('correlation').reset_index(drop=True)
s2.to_csv(OUT + 'set2_A_loose.csv', index=False)
print(f'  {len(s2):,} pairs')
del p, neg, pos, p_neg; gc.collect()

# ---------------------------------------------------------------- set1 and set2
print('\n[3/4] set1_inter_set2.csv')
inter = s2.merge(s1m.rename(columns={'miRNA': 'mirna_name', 'gene': 'gene_name'}),
                 on=['mirna_name', 'gene_name'], how='inner')
inter.sort_values('correlation').reset_index(drop=True).to_csv(OUT + 'set1_inter_set2.csv', index=False)
print(f'  {len(inter):,} pairs')

# ---------------------------------------------------------------- null: evaluable pairs in neither set
print('\n[4/4] null_relations_no_interaction.parquet')
gidx = {g: i for i, g in enumerate(G_arr)}
midx = {m: i for i, m in enumerate(M_arr)}
ex = np.zeros(len(rflat), dtype=bool)
for mi, ge in zip(s1m['miRNA'], s1m['gene']):
    ex[gidx[ge] * nM + midx[mi]] = True
for mi, ge in zip(s2['mirna_name'], s2['gene_name']):
    ex[gidx[ge] * nM + midx[mi]] = True

pos = np.flatnonzero(finite & ~ex)
rv = np.round(rflat[pos], 4).astype(np.float32)
del rflat, finite, ex; gc.collect()

schema = pa.schema([('mirna_name', pa.string()), ('gene_name', pa.string()), ('correlation', pa.float32())])
with pq.ParquetWriter(OUT + 'null_relations_no_interaction.parquet', schema) as writer:
    for s in range(0, len(pos), CHUNK):
        e = min(s + CHUNK, len(pos))
        writer.write_table(pa.table({'mirna_name': M_arr[pos[s:e] % nM].astype(str),
                                     'gene_name': G_arr[pos[s:e] // nM].astype(str),
                                     'correlation': rv[s:e]}, schema=schema))
print(f'  {len(pos):,} pairs')
