# ==============================================================================
# 03_compute_confusion_metrics.py - confusion matrix and error rates
# ==============================================================================
# Author : Congliang Jin
# Date   : 2026-10-23
#
# Description:
#   Truth = reference score >= c (c = 1, 2, 3). Prediction = top 10%, 20%, ..., 100% of
#   each tool's own ranked predictions (ties at the threshold included; pairs a tool
#   does not score count as not predicted). Writes TP/FP/FN/TN, FDP, Type I and Type II
#   error, sensitivity, specificity, Youden J, and the J-maximising cutoff of each tool.
#
# Input:
#   analysis/results/ref_score_table_evaluable.csv.gz
#
# Output (analysis/results/):
#   refscore_metrics.csv
#   refscore_best_cutoff.csv
#
# Prerequisite:
#   analysis/02_build_ref_score_table.py
#
# Memory:
#   large (reads the whole table)
#
# Run from the repository root:
#   python analysis/03_compute_confusion_metrics.py
# ==============================================================================
import os
import numpy as np
import pandas as pd


TABLE   = 'analysis/results/ref_score_table_evaluable.csv.gz'
OUT_DIR = 'analysis/results'
TOOLS = ['DIANA-microT', 'miRDB', 'mirDIP', 'miRmap', 'miRWalk', 'RNA22', 'TargetScan']

# False = a lower score is a stronger prediction (miRmap, RNA22 p-value, TargetScan context++)
HIGHER_BETTER = {'DIANA-microT': True, 'miRDB': True, 'mirDIP': True, 'miRWalk': True,
                 'miRmap': False,
                 'RNA22': False,
                 'TargetScan': False}
TRUTH_CUTOFFS = [1, 2, 3]
PCTS = list(range(10, 101, 10))


EXPECTED_REF_COUNTS = {0: 57_369_107, 1: 2_252_524, 2: 110_890, 3: 1_683}


def compute():
    print('Reading the table (a few minutes) ...', flush=True)
    ref_counts = np.zeros(4, dtype=np.int64)
    score_parts = {t: [] for t in TOOLS}
    ref_parts = {t: [] for t in TOOLS}
    reader = pd.read_csv(TABLE, usecols=['Ref_score', *TOOLS],
                         dtype={'Ref_score': 'int8', **{t: 'float64' for t in TOOLS}},
                         chunksize=2_000_000)
    for n_chunk, chunk in enumerate(reader, 1):
        ref = chunk['Ref_score'].to_numpy()
        ref_counts += np.bincount(ref, minlength=4)
        for t in TOOLS:
            v = chunk[t].to_numpy()
            has_score = ~np.isnan(v)
            score_parts[t].append(v[has_score]); ref_parts[t].append(ref[has_score])
        print(f'  {int(ref_counts.sum()):>12,} pairs read', flush=True)

    n_all = int(ref_counts.sum())
    print(f'\nPairs in the table: {n_all:,}')
    for s in (3, 2, 1, 0):
        note = ''
        if EXPECTED_REF_COUNTS is not None:
            note = '  (as expected)' if ref_counts[s] == EXPECTED_REF_COUNTS[s] else f'  (EXPECTED {EXPECTED_REF_COUNTS[s]:,})'
        print(f'  Ref score {s}: {int(ref_counts[s]):>12,}{note}')
    n_true = {c: int(ref_counts[c:].sum()) for c in TRUTH_CUTOFFS}

    rows, best_rows = [], []
    for t in TOOLS:
        v = np.concatenate(score_parts[t]); ref = np.concatenate(ref_parts[t])
        score_parts[t] = ref_parts[t] = None
        strength = v if HIGHER_BETTER[t] else -v
        order = np.argsort(-strength, kind='stable')
        strength, v, ref = strength[order], v[order], ref[order]
        n_pred_all = len(v)
        print(f'{t:<13} pairs with a score: {n_pred_all:>11,}', flush=True)

        k_distinct = np.append(np.flatnonzero(np.diff(strength) != 0) + 1, n_pred_all)

        for c in TRUTH_CUTOFFS:
            cum_tp = np.cumsum(ref >= c)
            P, NEG = n_true[c], n_all - n_true[c]

            for pct in PCTS:
                thr = np.quantile(strength, 1 - pct / 100)
                k = int(np.searchsorted(-strength, -thr, side='right'))
                tp = int(cum_tp[k - 1]) if k else 0
                fp, fn = k - tp, P - tp
                tn = NEG - fp
                rows.append(dict(tool=t, truth_cutoff=c, pct=pct,
                                 score_threshold=thr if HIGHER_BETTER[t] else -thr,
                                 n_predicted=k, TP=tp, FP=fp, FN=fn, TN=tn,
                                 FDP=fp / k if k else np.nan,
                                 type1=fp / NEG, type2=fn / P,
                                 sensitivity=tp / P, specificity=tn / NEG,
                                 youden=tp / P - fp / NEG))

            tp_k = cum_tp[k_distinct - 1]
            J = tp_k / P - (k_distinct - tp_k) / NEG
            b = int(np.argmax(J)); k = int(k_distinct[b]); tp = int(tp_k[b])
            best_rows.append(dict(tool=t, truth_cutoff=c,
                                  best_score_threshold=float(v[k - 1]),
                                  direction='score >= threshold' if HIGHER_BETTER[t] else 'score <= threshold',
                                  pct_of_predictions=100 * k / n_pred_all, n_predicted=k,
                                  TP=tp, FP=k - tp, FN=P - tp, TN=NEG - (k - tp),
                                  sensitivity=tp / P, specificity=1 - (k - tp) / NEG, youden=float(J[b])))
        del v, ref, strength, order
    return pd.DataFrame(rows), pd.DataFrame(best_rows)


def main():
    metrics, best = compute()

    best['FDP'] = best['FP'] / (best['TP'] + best['FP'])
    best['type1'] = best['FP'] / (best['FP'] + best['TN'])
    best['type2'] = best['FN'] / (best['TP'] + best['FN'])
    metrics.to_csv(os.path.join(OUT_DIR, 'refscore_metrics.csv'), index=False)
    best.to_csv(os.path.join(OUT_DIR, 'refscore_best_cutoff.csv'), index=False)

    show = best[['tool', 'truth_cutoff', 'best_score_threshold', 'direction', 'pct_of_predictions',
                 'sensitivity', 'specificity', 'youden', 'FDP', 'type1', 'type2']]
    with pd.option_context('display.width', 250, 'display.float_format', '{:.4g}'.format):
        print('\nBest cutoff per tool (maximum sensitivity + specificity)')
        print(show.to_string(index=False))
    print(f'\nSaved refscore_metrics.csv and refscore_best_cutoff.csv in {OUT_DIR}')


if __name__ == '__main__':
    main()
