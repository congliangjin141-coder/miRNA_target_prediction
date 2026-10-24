# ==============================================================================
# verify_ref_table.py - independent check of the reference table
# ==============================================================================
# Author : Congliang Jin
# Date   : 2026-10-23
#
# Description:
#   Every column of the table is recomputed from the input files and compared pair by pair:
#     1. structure   one row per candidate pair; 2,036 x 29,339 = 59,734,204 rows
#     2. iPSC        Y exactly for the pairs of set2_A_loose.csv (r < -0.2, p < 0.05, FDR < 0.05)
#     3. miRTarBase  Y exactly for the pairs of set1_in_matrix.csv with source miRTarBase or both
#     4. TarBase     Y exactly for the pairs of set1_in_matrix.csv with source TarBase or both
#     5. Ref_score   the number of Y among the three references
#     6. score 0     exactly the pairs of null_relations_no_interaction.parquet
#     7. tools       the strongest score of each pair in each tool file
#   The checks use their own code (pandas group-by on pair keys), not the functions of
#   02_build_ref_score_table.py.
#
# Input:
#   analysis/results/ref_score_table_evaluable.csv.gz
#   data/reference/ and data/predictions/
#
# Output (analysis/results/):
#   verify_ref_table_report.txt
#
# Prerequisite:
#   analysis/02_build_ref_score_table.py
#
# Memory:
#   about 64 GB; 1-3 hours
#
# Run from the repository root:
#   python analysis/checks/verify_ref_table.py
# ==============================================================================
import numpy as np
import pandas as pd

REF = 'data/reference'
PRED = 'data/predictions'
TABLE = 'analysis/results/ref_score_table_evaluable.csv.gz'
REPORT = 'analysis/results/verify_ref_table_report.txt'

SET1 = f'{REF}/set1_in_matrix.csv'                      # miRNA, gene, source
SET2 = f'{REF}/set2_A_loose.csv'                        # mirna_name, gene_name, correlation, pvalue, FDR
NULL = f'{REF}/null_relations_no_interaction.parquet'   # mirna_name, gene_name
TOOLS = {  # column in the table: (file, how several rows of one pair are combined)
    'TargetScan':   ('TargetScan v8.0.csv', 'min'),
    'RNA22':        ('rna22.csv', 'min'),
    'DIANA-microT': ('diana_microT.csv', 'max'),
    'miRmap':       ('mirmap.csv', 'min'),
    'miRDB':        ('miRDB v6.0.csv', 'max'),
    'mirDIP':       ('mirDIP v5.2.csv', 'max'),
    'miRWalk':      ('mirwalk.csv', 'max'),
}
EXPECTED_SHAPE = (2036, 29339)
CHUNK = 5_000_000

lines, results = [], []


def report(text=''):
    print(text, flush=True)
    lines.append(text)


def check(label, ok, detail=''):
    results.append(ok)
    report(f'  [{"OK  " if ok else "FAIL"}] {label}' + (f'   ({detail})' if detail else ''))


# ----------------------------------------------------------------- 1. read the table
report('1. Reading the table ...')
mir_codes, gene_codes, cols = [], [], {c: [] for c in ['miRTarBase', 'TarBase', 'iPSC', 'Ref_score', *TOOLS]}
mir_names, gene_names = {}, {}
bad_yn = {'miRTarBase': 0, 'TarBase': 0, 'iPSC': 0}
for chunk in pd.read_csv(TABLE, dtype={'miRNA': str, 'Gene': str, 'miRTarBase': str, 'TarBase': str,
                                       'iPSC': str, 'Ref_score': 'int8',
                                       **{t: 'float64' for t in TOOLS}}, chunksize=CHUNK):
    for names, store, col in ((mir_names, mir_codes, 'miRNA'), (gene_names, gene_codes, 'Gene')):
        new = pd.unique(chunk[col][~chunk[col].isin(names.keys())])
        for n in new:
            names[n] = len(names)
        store.append(chunk[col].map(names).to_numpy(np.int64))
    for c in cols:
        if c in ('miRTarBase', 'TarBase', 'iPSC'):
            bad_yn[c] += int((~chunk[c].isin(['Y', 'N'])).sum())
            cols[c].append((chunk[c] == 'Y').to_numpy())
        else:
            cols[c].append(chunk[c].to_numpy())
    report(f'    {sum(len(x) for x in mir_codes):>12,} rows')

mir_codes, gene_codes = np.concatenate(mir_codes), np.concatenate(gene_codes)
cols = {c: np.concatenate(v) for c, v in cols.items()}
n_rows, n_mir, n_gene = len(mir_codes), len(mir_names), len(gene_names)
key = mir_codes * n_gene + gene_codes                      # one integer per pair
report(f'\n1. Structure')
check('number of miRNAs x genes', (n_mir, n_gene) == EXPECTED_SHAPE, f'{n_mir:,} x {n_gene:,}')
check('number of rows = miRNAs x genes', n_rows == n_mir * n_gene, f'{n_rows:,} rows')
check('every pair appears exactly once', np.bincount(key, minlength=n_mir * n_gene).max() == 1)
for c in ('miRTarBase', 'TarBase', 'iPSC'):
    check(f'{c} column holds only Y / N', bad_yn[c] == 0)
order = np.argsort(key)                                    # table values in pair-key order
flag = {c: cols[c][order] for c in ('miRTarBase', 'TarBase', 'iPSC')}
ref_score = cols['Ref_score'][order]
tool_col = {t: cols[t][order] for t in TOOLS}
del cols, mir_codes, gene_codes, key, order


def keys_of(mirna, gene):
    """Pair key of each (miRNA, gene); -1 when the miRNA or the gene is not in the table."""
    m = pd.Series(mirna).map(mir_names)
    g = pd.Series(gene).map(gene_names)
    inside = m.notna().to_numpy() & g.notna().to_numpy()
    k = np.full(len(m), -1, dtype=np.int64)
    k[inside] = m[inside].to_numpy(np.int64) * n_gene + g[inside].to_numpy(np.int64)
    return k


def as_mask(keys):
    mask = np.zeros(n_mir * n_gene, dtype=bool)
    mask[keys] = True
    return mask


# ----------------------------------------------------------------- 2. iPSC
report('\n2. iPSC (set2_A_loose.csv)')
s2 = pd.read_csv(SET2, dtype={'mirna_name': str, 'gene_name': str})
s2['mirna_name'], s2['gene_name'] = s2['mirna_name'].str.strip(), s2['gene_name'].str.strip()
check('every pair meets r < -0.2, p < 0.05, FDR < 0.05',
      bool(((s2['correlation'] < -0.2) & (s2['pvalue'] < 0.05) & (s2['FDR'] < 0.05)).all()), f'{len(s2):,} pairs')
k2 = keys_of(s2['mirna_name'], s2['gene_name'])
check('every iPSC pair is a candidate pair', bool((k2 >= 0).all()))
check('table iPSC = Y exactly for these pairs', np.array_equal(as_mask(k2[k2 >= 0]), flag['iPSC']),
      f'table Y {int(flag["iPSC"].sum()):,}')
del s2, k2

# ----------------------------------------------------------------- 3-4. databases
report('\n3-4. miRTarBase / TarBase (set1_in_matrix.csv)')
s1 = pd.read_csv(SET1, dtype=str)
s1['miRNA'], s1['gene'] = s1['miRNA'].str.strip(), s1['gene'].str.strip()
src = s1['source'].str.strip().str.lower()
report(f'    source values: {dict(src.value_counts())}')
k1 = keys_of(s1['miRNA'], s1['gene'])
report(f'    database rows {len(s1):,}; candidate pairs {int((k1 >= 0).sum()):,}; '
       f'outside the candidate set {int((k1 < 0).sum()):,}')
for name, wanted in (('miRTarBase', {'mirtarbase', 'both'}), ('TarBase', {'tarbase', 'both'})):
    sel = src.isin(wanted).to_numpy() & (k1 >= 0)
    check(f'table {name} = Y exactly for its candidate pairs', np.array_equal(as_mask(k1[sel]), flag[name]),
          f'table Y {int(flag[name].sum()):,}')
del s1, k1, src

# ----------------------------------------------------------------- 5-6. reference score and score 0
report('\n5-6. Reference score')
check('Ref_score = number of Y', np.array_equal(
    ref_score, flag['miRTarBase'].astype(np.int8) + flag['TarBase'] + flag['iPSC']))
report('    pairs by Ref_score: ' + ', '.join(f'{s}: {int((ref_score == s).sum()):,}' for s in (3, 2, 1, 0)))
import pyarrow.parquet as pq
nul = pq.read_table(NULL, columns=['mirna_name', 'gene_name']).to_pandas()
kn = keys_of(nul['mirna_name'].astype(str), nul['gene_name'].astype(str))
check('every null pair is a candidate pair', bool((kn >= 0).all()), f'{len(nul):,} null pairs')
check('Ref_score = 0 exactly for the null pairs', np.array_equal(as_mask(kn[kn >= 0]), ref_score == 0))
del nul, kn, flag, ref_score

# ----------------------------------------------------------------- 7. tools
report('\n7. Tool scores')
for name, (fname, how) in TOOLS.items():
    parts = []
    for chunk in pd.read_csv(f'{PRED}/{fname}', usecols=['mirna_name', 'gene_name', 'score'], dtype=str,
                             chunksize=CHUNK):
        chunk['score'] = pd.to_numeric(chunk['score'], errors='coerce')
        chunk['k'] = keys_of(chunk['mirna_name'], chunk['gene_name'])
        chunk = chunk[(chunk['k'] >= 0) & chunk['score'].notna()]
        parts.append(chunk.groupby('k')['score'].agg(how))
        if len(parts) >= 20:
            parts = [pd.concat(parts).groupby(level=0).agg(how)]
    best = pd.concat(parts).groupby(level=0).agg(how)
    expected = np.full(n_mir * n_gene, np.nan)
    expected[best.index.to_numpy()] = best.to_numpy()
    observed = tool_col[name]
    same_mask = np.array_equal(np.isnan(expected), np.isnan(observed))
    both = ~np.isnan(expected) & ~np.isnan(observed)
    max_diff = float(np.max(np.abs(expected[both] - observed[both]))) if both.any() else 0.0
    check(f'{name:<13} same scored pairs', same_mask,
          f'file {int((~np.isnan(expected)).sum()):,} / table {int((~np.isnan(observed)).sum()):,}')
    check(f'{name:<13} same strongest score', max_diff <= 1e-12 * max(1.0, float(np.nanmax(np.abs(observed)))),
          f'largest difference {max_diff:.3g}')
    del parts, best, expected

report('\nRESULT: ' + ('the table matches its input files in every check' if all(results)
                       else f'{results.count(False)} check(s) FAILED, see above'))
with open(REPORT, 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print(f'Report saved to {REPORT}')
