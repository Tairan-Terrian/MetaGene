"""Recheck evidence counts without changing predictions, thresholds or candidate ranks."""
import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parents[1]


def read_rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def review_evidence(root=ROOT):
    directory = root / 'results/candidates'
    core = read_rows(directory / 'core_candidates.csv')
    top500 = read_rows(directory / 'top500.csv')
    top30 = read_rows(directory / 'top30.csv')
    representatives = read_rows(directory / 'representative3.csv')
    top_ranks = {row['pair_key']: int(row['top500_rank']) for row in top500}
    short_pairs = {row['pair_key'] for row in top30}
    cutoff = json.loads((directory / 'ranking_policy.json').read_text())['single_cutoff']
    recomputed = {}
    followup = []
    for row in core:
        associated_cells, qualifying_arms = [], []
        complete = 0
        for cell in ['HCT116', 'SW480']:
            associated = False
            for arm in ['a', 'b']:
                for dose in ['HD', 'MD']:
                    prefix = f'gene_{arm}_GSE299428_{cell}_{dose}_vs_PAR_'
                    effect, fdr = numeric(row[prefix + 'log2fc']), numeric(row[prefix + 'padj'])
                    complete += int(math.isfinite(effect) and math.isfinite(fdr))
                    if fdr < 0.05 and abs(effect) >= 1:
                        associated = True
                        qualifying_arms.append((row['gene_' + arm], cell, dose, effect, fdr))
            if associated:
                associated_cells.append(cell)
        if len(associated_cells) != int(row['target_resistance_association_count']):
            raise ValueError(f'Resistance-expression count mismatch: {row["pair_key"]}')
        single_values = [numeric(row[f'gene_{arm}_colo1_single_control_scaled']) for arm in ['a', 'b']]
        observed_singles = sum(math.isfinite(x) for x in single_values)
        expected_strong = sum(x <= cutoff for x in single_values) if observed_singles == 2 else None
        actual_strong = numeric(row['strong_single_count'])
        if (expected_strong is None and math.isfinite(actual_strong)) or (expected_strong is not None and actual_strong != expected_strong):
            raise ValueError(f'Single-effect count mismatch: {row["pair_key"]}')
        if int(row['external_measurement_count']) == 0 and row['external_pair_state'] != 'not_covered_unknown':
            raise ValueError(f'Unmeasured external pair mislabeled: {row["pair_key"]}')
        recomputed[row['pair_key']] = {'associated_cells': len(associated_cells), 'complete_expression_comparisons': complete}
        if associated_cells and row['pair_key'] in top_ranks:
            followup.append({
                'top500_rank': top_ranks[row['pair_key']], 'pair_key': row['pair_key'],
                'metabolic_gene': row['metabolic_gene'], 'histone_enzyme_gene': row['chromatin_gene'],
                'pred_gi': row['pred_gi'], 'pred_joint': row['pred_joint'],
                'observed_HT29_single_count': row['observed_HT29_single_count'],
                'missing_HT29_single_genes': ';'.join(row['gene_' + a] for a in ['a', 'b'] if not math.isfinite(numeric(row['single_' + a]))),
                'associated_cell_backgrounds': ';'.join(associated_cells),
                'qualifying_gene_evidence': ';'.join(f'{g}:{c}_{d}:log2FC={e:.6g}:FDR={p:.6g}' for g, c, d, e, p in qualifying_arms),
                'in_top30': row['pair_key'] in short_pairs,
                'shortlist_status': 'excluded_missing_HT29_single_measurements' if int(row['observed_HT29_single_count']) < 2 else 'eligible_for_evidence_ranking',
                'interpretation': 'Gene-level expression association only; not pair-specific efficacy or resistance reversal; no change to the frozen shortlist.',
            })
    followup.sort(key=lambda r: r['top500_rank'])
    curves = read_rows(root / 'results/active_learning/curves.csv')
    result = {
        'reviewed_at': datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='minutes'),
        'scope': 'Recomputed evidence counts from saved results; no new fitting, threshold changes or reranking.',
        'core_pairs': len(core), 'top500_pairs': len(top500), 'top30_pairs': len(top30),
        'representative_pairs': len(representatives),
        'top500_with_target_expression_association': len(followup),
        'associated_pairs_excluded_by_missing_HT29_singles': sum(int(r['observed_HT29_single_count']) < 2 for r in followup),
        'associated_metabolic_genes': sorted({r['metabolic_gene'] for r in followup}),
        'top30_with_target_expression_association': sum(recomputed[r['pair_key']]['associated_cells'] > 0 for r in top30),
        'representatives_with_target_expression_association': sum(recomputed[r['pair_key']]['associated_cells'] > 0 for r in representatives),
        'top30_expression_comparison_coverage': dict(Counter(recomputed[r['pair_key']]['complete_expression_comparisons'] for r in top30)),
        'top30_external_measured_pairs': sum(int(r['external_measurement_count']) > 0 for r in top30),
        'top30_external_hit_rate': (sum(int(r['external_hit_count']) > 0 for r in top30) / sum(int(r['external_measurement_count']) > 0 for r in top30)) if any(int(r['external_measurement_count']) > 0 for r in top30) else None,
        'top30_external_hit_rate_reason': 'Computed only among independently measured candidate pairs; undefined with no covered candidate. Absence is not a negative result.',
        'top30_without_strong_single_effect': sum(numeric(r['strong_single_count']) == 0 for r in top30),
        'active_learning_core_measured_pool': int(curves[0]['core_pool_n']),
        'active_learning_core_hit_k': int(curves[0]['core_hit_k']),
        'delivery_assessment': 'Computational-stage report and exploratory candidate handoff; insufficient evidence for a finalized set of reliably predicted resistance-reversal combinations.',
    }
    output = root / 'results/handoff'
    output.mkdir(parents=True, exist_ok=True)
    columns = ['top500_rank', 'pair_key', 'metabolic_gene', 'histone_enzyme_gene', 'pred_gi', 'pred_joint',
               'observed_HT29_single_count', 'missing_HT29_single_genes', 'associated_cell_backgrounds',
               'qualifying_gene_evidence', 'in_top30', 'shortlist_status', 'interpretation']
    with (output / 'resistance_expression_followup.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(followup)
    (output / 'evidence_review.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    root=parser.parse_args().root
    current=root/'results/sl/delivery_review.json'
    review=json.loads(current.read_text()) if current.exists() else review_evidence(root)
    print(json.dumps(review, ensure_ascii=False, indent=2))
