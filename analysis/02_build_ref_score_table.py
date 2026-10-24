# ==============================================================================
# 02_build_ref_score_table.py - combined reference table
# ==============================================================================
# Author : Congliang Jin
# Date   : 2026-10-23
#
# Description:
#   One row per candidate pair (2,036 miRNAs x 29,339 genes = 59,734,204 pairs) with
#   Y/N for miRTarBase, TarBase and iPSC, the reference score (0-3) and the score of
#   each of the seven tools. Several rows of one pair in a tool file are combined into
#   the strongest score (minimum for TargetScan, RNA22, miRmap; maximum otherwise).
#
# Input:
#   data/reference/  set1_in_matrix.csv, set2_A_loose.csv, null_relations_no_interaction.parquet
#   data/predictions/  the seven tool files (see data/README.md)
#
# Output (analysis/results/):
#   ref_score_table_evaluable.csv.gz
#   ref_score_summary_evaluable.csv
#
# Prerequisite:
#   analysis/01_generate_input_sets.py
#
# Memory:
#   about 64 GB
#
# Run from the repository root:
#   python analysis/02_build_ref_score_table.py
# ==============================================================================
import gzip
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(".")
REF = BASE / "data" / "reference"
PRED = BASE / "data" / "predictions"
OUT_DIR = BASE / "analysis" / "results"
OUT_TABLE = "ref_score_table_evaluable.csv.gz"
OUT_SUMMARY = "ref_score_summary_evaluable.csv"


SET1 = dict(path=REF / "set1_in_matrix.csv", sep=",", mirna="miRNA", gene="gene", source="source")


IPSC_SIG = dict(path=REF / "set2_A_loose.csv", sep=",", mirna="mirna_name", gene="gene_name",
                corr="correlation", pvalue="pvalue", fdr="FDR")
# iPSC = Y: Pearson r < -0.2, p < 0.05 and BH FDR < 0.05
R_CUT, P_CUT, FDR_CUT = -0.2, 0.05, 0.05


NULL = dict(path=REF / "null_relations_no_interaction.parquet", sep=",",
            mirna="mirna_name", gene="gene_name")


# several rows of one pair are combined into its strongest score: min = lower is stronger
METHODS = {
    "TargetScan":   dict(path=PRED / "TargetScan v8.0.csv", sep=",", mirna="mirna_name", gene="gene_name", score="score", agg="min"),
    "RNA22":        dict(path=PRED / "rna22.csv",           sep=",", mirna="mirna_name", gene="gene_name", score="score", agg="min"),
    "DIANA-microT": dict(path=PRED / "diana_microT.csv",    sep=",", mirna="mirna_name", gene="gene_name", score="score", agg="max"),
    "miRmap":       dict(path=PRED / "mirmap.csv",          sep=",", mirna="mirna_name", gene="gene_name", score="score", agg="min"),
    "miRDB":        dict(path=PRED / "miRDB v6.0.csv",      sep=",", mirna="mirna_name", gene="gene_name", score="score", agg="max"),
    "mirDIP":       dict(path=PRED / "mirDIP v5.2.csv",     sep=",", mirna="mirna_name", gene="gene_name", score="score", agg="max"),
    "miRWalk":      dict(path=PRED / "mirwalk.csv",         sep=",", mirna="mirna_name", gene="gene_name", score="score", agg="max"),
}
EXPECTED_SHAPE = (2036, 29339)

CHUNK_ROWS = 5_000_000


EXPECTED = {
    (1, 1, 1): 1_683, (1, 1, 0): None, (1, 0, 1): 6_013, (0, 1, 1): 24_041,
    (1, 0, 0): None, (0, 1, 0): None, (0, 0, 1): 225_139, (0, 0, 0): 57_369_107,
}


REF_NAMES = ["miRTarBase", "TarBase", "iPSC"]
YN = np.array(["N", "Y"])


def read_pairs(cfg, extra=None):
    extra = {k: v for k, v in (extra or {}).items() if k is not None}
    df = pd.read_csv(cfg["path"], sep=cfg["sep"], usecols=[cfg["mirna"], cfg["gene"], *extra],
                     dtype={cfg["mirna"]: str, cfg["gene"]: str})
    df = df.rename(columns={cfg["mirna"]: "miRNA", cfg["gene"]: "Gene", **extra})
    df["miRNA"] = df["miRNA"].str.strip()
    df["Gene"] = df["Gene"].str.strip()
    return df.dropna(subset=["miRNA", "Gene"]).reset_index(drop=True)


def read_null(cfg):
    names = [cfg["mirna"], cfg["gene"]]
    path = Path(cfg["path"])
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        df = pq.read_table(path, columns=names, read_dictionary=names).to_pandas()
    else:
        df = pd.read_csv(path, sep=cfg["sep"], usecols=names, dtype={n: "category" for n in names})
    df = df.rename(columns={cfg["mirna"]: "miRNA", cfg["gene"]: "Gene"})
    df = df.dropna(subset=["miRNA", "Gene"]).reset_index(drop=True)
    for c in ["miRNA", "Gene"]:
        df[c] = df[c].astype("category").cat.remove_unused_categories()
    return df


def source_flags(source):
    s = source.fillna("").str.lower()
    both = s.str.contains("both")
    in_mtb = s.str.contains("mirtarbase") | both
    in_tb = s.str.replace("mirtarbase", "", regex=False).str.contains("tarbase") | both
    if (~(in_mtb | in_tb)).any():
        raise ValueError(f"Unrecognised `source` values: {sorted(source[~(in_mtb | in_tb)].unique())[:10]}")
    return in_mtb.to_numpy(), in_tb.to_numpy()


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    set1 = read_pairs(SET1, extra={SET1["source"]: "source"})
    print("set1 `source` values:\n" + set1["source"].value_counts(dropna=False).to_string())
    sig = read_pairs(IPSC_SIG, extra={IPSC_SIG["corr"]: "corr", IPSC_SIG["pvalue"]: "p", IPSC_SIG["fdr"]: "fdr"})
    null = read_null(NULL)
    print(f"\nset1 rows: {len(set1):,}   iPSC = Y rows: {len(sig):,}   null rows: {len(null):,}")

    mirs = pd.Index(sorted(null["miRNA"].cat.categories))
    genes = pd.Index(sorted(null["Gene"].cat.categories))
    shape = (len(mirs), len(genes))
    n_cells = shape[0] * shape[1]
    print(f"Evaluable set: {shape[0]:,} miRNAs x {shape[1]:,} genes = {n_cells:,} pairs")
    if EXPECTED_SHAPE is not None and shape != tuple(EXPECTED_SHAPE):
        raise ValueError(f"Expected {EXPECTED_SHAPE[0]:,} miRNAs x {EXPECTED_SHAPE[1]:,} genes.")

    def locate(df):
        i, j = mirs.get_indexer(df["miRNA"]), genes.get_indexer(df["Gene"])
        return i, j, (i >= 0) & (j >= 0)

    i = mirs.get_indexer(null["miRNA"].cat.categories)[null["miRNA"].cat.codes.to_numpy()]
    j = genes.get_indexer(null["Gene"].cat.categories)[null["Gene"].cat.codes.to_numpy()]
    in_null = np.zeros(shape, dtype=bool)
    in_null[i, j] = True
    if int(in_null.sum()) != len(null):
        raise ValueError("The null file has duplicated pairs.")
    del null, i, j

    if sig.duplicated(["miRNA", "Gene"]).any():
        raise ValueError("set2 file has duplicated pairs; expected one row per pair.")
    bad = ~(sig["corr"] < R_CUT)
    if "p" in sig:
        bad |= ~(sig["p"] < P_CUT)
    if "fdr" in sig:
        bad |= ~(sig["fdr"] < FDR_CUT)
    if bad.any():
        raise ValueError(f"{int(bad.sum()):,} pairs of the set2 file do not meet r < {R_CUT}, p < {P_CUT}, FDR < {FDR_CUT}.")
    i, j, inside = locate(sig)
    if not inside.all():
        raise ValueError(f"{int((~inside).sum()):,} iPSC = Y pairs involve a miRNA or gene that is not in the evaluable set.")
    flags = {"iPSC": np.zeros(shape, dtype=bool)}
    flags["iPSC"][i, j] = True
    print(f"iPSC = Y: {len(sig):,} pairs, all meet r < {R_CUT}, p < {P_CUT}, FDR < {FDR_CUT}")

    in_mtb, in_tb = source_flags(set1["source"])
    i, j, inside = locate(set1)
    print(f"Database pairs: {len(set1):,}   evaluable (kept): {int(inside.sum()):,}   not evaluable (left out): {int((~inside).sum()):,}")
    for label, sel in [("miRTarBase only", in_mtb & ~in_tb), ("TarBase only", in_tb & ~in_mtb), ("both databases", in_mtb & in_tb)]:
        print(f"  left out, {label:<16}: {int((sel & ~inside).sum()):,}")
    for name, keep in [("miRTarBase", in_mtb), ("TarBase", in_tb)]:
        flags[name] = np.zeros(shape, dtype=bool)
        flags[name][i[keep & inside], j[keep & inside]] = True
        print(f"{name:<11} evaluable pairs: {int(flags[name].sum()):,}")
    del set1, sig

    in_ref = flags["miRTarBase"] | flags["TarBase"] | flags["iPSC"]
    n_overlap = int((in_null & in_ref).sum())
    if n_overlap:
        raise ValueError(f"{n_overlap:,} null pairs are also in a reference -> the null file is not the score-0 set.")
    n_missing = int((~(in_null | in_ref)).sum())
    if n_missing:
        raise ValueError(f"{n_missing:,} evaluable pairs are neither in the null file nor in a reference.")
    print(f"Coverage check passed: null ({int(in_null.sum()):,}) + references ({int(in_ref.sum()):,}) = {n_cells:,} pairs, no overlap.")
    del in_null, in_ref

    score = sum(flags[n].astype(np.int8) for n in REF_NAMES)

    rows = []
    for a, b, c in product([1, 0], repeat=3):
        mask = (flags["miRTarBase"] == a) & (flags["TarBase"] == b) & (flags["iPSC"] == c)
        row = {"Ref_score": a + b + c, "miRTarBase": YN[a], "TarBase": YN[b], "iPSC": YN[c],
               "n_pairs": int(mask.sum()),
               "n_miRNA": int(mask.any(axis=1).sum()),
               "n_Gene": int(mask.any(axis=0).sum())}
        if EXPECTED is not None:
            want = EXPECTED[(a, b, c)]
            row["n_pairs_expected"] = "" if want is None else want
            row["match"] = "n/a" if want is None else ("yes" if row["n_pairs"] == want else "NO")
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values("Ref_score", ascending=False, kind="stable")
    summary.to_csv(OUT_DIR / OUT_SUMMARY, index=False)
    print("\nVenn-region summary")
    print(summary.to_string(index=False))
    if EXPECTED is not None and (summary["match"] == "NO").any():
        print("WARNING: some regions differ from the expected counts.")

    scores = {}
    for name, cfg in METHODS.items():
        best = np.fmin if cfg["agg"] == "min" else np.fmax
        M = np.full(shape, np.nan)
        flat_M = M.reshape(-1)
        n_rows = n_outside = n_bad = 0

        reader = pd.read_csv(cfg["path"], sep=cfg["sep"], usecols=[cfg["mirna"], cfg["gene"], cfg["score"]],
                             dtype=str, chunksize=CHUNK_ROWS)
        for chunk in reader:
            chunk = chunk.rename(columns={cfg["mirna"]: "miRNA", cfg["gene"]: "Gene", cfg["score"]: "score"})
            value = pd.to_numeric(chunk["score"], errors="coerce").to_numpy(dtype=float)
            i, j, inside = locate(chunk)
            n_rows += len(chunk)
            n_bad += int(np.isnan(value).sum())
            n_outside += int((~inside).sum())
            inside &= ~np.isnan(value)
            flat = i[inside].astype(np.int64) * shape[1] + j[inside]
            per_pair = pd.Series(value[inside]).groupby(flat).agg(cfg["agg"])
            idx = per_pair.index.to_numpy()
            flat_M[idx] = best(flat_M[idx], per_pair.to_numpy())
        scores[name] = M
        print(f"{name:<13} rows read: {n_rows:>13,}   rows not evaluable (dropped): {n_outside:>12,}   "
              f"rows without a numeric score (dropped): {n_bad:>9,}   "
              f"pairs with a score: {int((~np.isnan(M)).sum()):>11,}", flush=True)

    cols = ["miRNA", "Gene", *REF_NAMES, "Ref_score", *METHODS]
    path = OUT_DIR / OUT_TABLE
    fh = gzip.open(path, "wt", newline="", compresslevel=4) if path.suffix == ".gz" else open(path, "w", newline="")
    gene_values = genes.to_numpy()
    with fh:
        for r, mir in enumerate(mirs):
            block = {"miRNA": mir, "Gene": gene_values}
            block.update({n: YN[flags[n][r].astype(int)] for n in REF_NAMES})
            block["Ref_score"] = score[r]
            block.update({n: scores[n][r] for n in METHODS})
            pd.DataFrame(block)[cols].to_csv(fh, header=(r == 0), index=False)
    print(f"\nRows written to {path}: {n_cells:,}")


if __name__ == "__main__":
    main()
