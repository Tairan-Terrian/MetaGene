"""Independent DepMap, TCGA, single-cell and STRING covariates for all model genes.

No double-perturbation labels are used. Statistical normalization for training
is fitted later on each training gene set. Missing covariates stay NaN and have
explicit indicator columns. Patient summaries weight patients equally.
"""
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/evidence'
SRC = OUT / 'requested_sources'


def depmap(genes):
    model = pd.read_csv(SRC / 'DepMap_Model.csv', index_col='ModelID')
    header = pd.read_csv(SRC / 'DepMap_CRISPRGeneEffect.csv', nrows=0).columns
    mapping = {column: re.sub(r' \(\d+\)$', '', column) for column in header[1:]}
    wanted = [column for column, gene in mapping.items() if gene in genes]
    raw = pd.read_csv(SRC / 'DepMap_CRISPRGeneEffect.csv', usecols=[header[0]] + wanted, index_col=0)
    raw = raw.rename(columns=mapping).T.groupby(level=0).mean().T
    annotated = model.reindex(raw.index)
    colorectal = annotated.OncotreePrimaryDisease.eq('Colorectal Adenocarcinoma')
    evidence = pd.DataFrame(index=genes)
    evidence['DepMap_pan_cancer_median'] = raw.median()
    evidence['DepMap_CRC_median'] = raw.loc[colorectal].median()
    evidence['DepMap_CRC_sd'] = raw.loc[colorectal].std()
    evidence['DepMap_CRC_n'] = raw.loc[colorectal].count()
    for cell in ['HT29', 'HCT116', 'SW480', 'K562', 'JURKAT']:
        ids = annotated.index[annotated.StrippedCellLineName.eq(cell)]
        evidence[f'DepMap_{cell}_gene_effect'] = raw.loc[ids].median() if len(ids) else np.nan
    keep = colorectal | annotated.StrippedCellLineName.isin(['HT29', 'HCT116', 'SW480', 'K562', 'JURKAT'])
    raw.loc[keep].to_csv(OUT / 'depmap_selected_cell_gene_effects.csv')
    annotated.loc[keep].to_csv(OUT / 'depmap_selected_models.csv')
    evidence.to_csv(OUT / 'depmap_gene_evidence.csv', index_label='gene')
    metadata = {'release': 'DepMap 24Q4 Public', 'doi': '10.25452/figshare.plus.27993248.v1',
                'total_screened_models': len(raw), 'CRC_screened_models': int(colorectal.sum()),
                'covered_model_genes': len(set(raw.columns) & set(genes)),
                'scope': 'Independent CRISPR single-gene dependency, parental cell models; not acquired-resistant HCT116-R/SW480-R and not combination effects.'}
    return evidence.drop(columns='DepMap_CRC_n'), metadata


def tcga(genes):
    raw = pd.read_csv(SRC / 'TCGA_COAD_expression.tsv.gz', sep='\t', index_col=0)
    raw.index = raw.index.astype(str)
    raw = raw.groupby(level=0).mean().reindex(genes)
    types = pd.Series({sample: sample.split('-')[3][:2] if len(sample.split('-')) >= 4 else '' for sample in raw.columns})
    rows = []
    evidence = pd.DataFrame(index=genes)
    counts = {}
    for name, code in [('tumor', '01'), ('normal', '11')]:
        columns = types.index[types.eq(code)]
        patient = raw[columns].T.groupby(lambda sample: '-'.join(sample.split('-')[:3])).median()
        counts[name] = {'samples': len(columns), 'patients': len(patient)}
        evidence[f'TCGA_{name}_median_log2_RSEM'] = patient.median()
        evidence[f'TCGA_{name}_fraction_above_log2_1'] = (patient.gt(1).sum() / patient.notna().sum()).replace([np.inf, -np.inf], np.nan)
        for sample in columns:
            rows.append({'sample': sample, 'patient': '-'.join(sample.split('-')[:3]), 'type': name})
    evidence['TCGA_tumor_minus_normal_median'] = evidence.TCGA_tumor_median_log2_RSEM - evidence.TCGA_normal_median_log2_RSEM
    raw.to_csv(OUT / 'tcga_model_gene_expression.csv', index_label='gene')
    pd.DataFrame(rows).to_csv(OUT / 'tcga_sample_manifest.csv', index=False)
    evidence.to_csv(OUT / 'tcga_gene_evidence.csv', index_label='gene')
    return evidence, {'cohort': 'TCGA-COAD', 'unit': 'log2(RSEM normalized count + 1)',
                      'counts': counts, 'covered_model_genes': int(raw.notna().any(axis=1).sum()),
                      'scope': 'Tumor/normal expression context; treatment response and oxaliplatin resistance are not inferred.'}


def single_cell(genes):
    cache = OUT / 'single_cell_model_counts.npz'
    annotation = pd.read_csv(SRC / 'GSE132465_cell_annotation.txt.gz', sep='\t', index_col='Index')
    if cache.exists():
        z = np.load(cache)
        assert z['genes'].tolist() == list(genes)
        counts, cells, total = z['counts'], z['cells'], z['library_size']
    else:
        index = {gene: i for i, gene in enumerate(genes)}
        with gzip.open(SRC / 'GSE132465_counts.txt.gz', 'rt') as stream:
            header = stream.readline().strip().split('\t')
            cells = np.asarray(header[1:], dtype=str)
            # Deposited matrix may omit the empty row-name field in its header.
            if len(cells) != len(annotation) and len(header) == len(annotation):
                cells = np.asarray(header, dtype=str)
            total = np.zeros(len(cells), dtype=np.float64)
            counts = np.zeros((len(genes), len(cells)), dtype=np.float32)
            present = np.zeros(len(genes), dtype=bool)
            for number, line in enumerate(stream, start=1):
                gene, values = line.rstrip().split('\t', 1)
                values = np.fromstring(values, sep='\t', dtype=np.float32)
                if len(values) != len(cells):
                    raise ValueError(f'Unexpected cell count in matrix row {number}')
                total += values
                if gene in index:
                    counts[index[gene]] += values
                    present[index[gene]] = True
                if number % 5000 == 0:
                    print(json.dumps({'single_cell_genes_read': number}), flush=True)
            counts[~present] = np.nan
        np.savez_compressed(cache, genes=np.asarray(genes, dtype=str), cells=cells, counts=counts, library_size=total)
    annotation = annotation.reindex(cells)
    assert annotation.Patient.notna().all() and (total > 0).all()
    norm = np.log1p(counts / total[None, :] * 10000).astype(np.float32)
    groups = {
        'tumor_epithelial': annotation.Class.eq('Tumor') & annotation.Cell_type.eq('Epithelial cells'),
        'normal_epithelial': annotation.Class.eq('Normal') & annotation.Cell_type.eq('Epithelial cells'),
        'tumor_nonepithelial': annotation.Class.eq('Tumor') & annotation.Cell_type.ne('Epithelial cells'),
    }
    evidence = pd.DataFrame(index=genes)
    metadata, patient_frames = {}, []
    for name, mask in groups.items():
        patients = sorted(annotation.loc[mask, 'Patient'].unique())
        means, fractions = [], []
        for patient in patients:
            select = (mask & annotation.Patient.eq(patient)).to_numpy()
            mean = np.nanmean(norm[:, select], axis=1)
            fraction = np.where(np.isfinite(counts[:, select]).any(axis=1), np.mean(counts[:, select] > 0, axis=1), np.nan)
            means.append(mean)
            fractions.append(fraction)
            patient_frames.append(pd.DataFrame({'gene': genes, 'patient': patient, 'compartment': name,
                                                'n_cells': int(select.sum()), 'mean_log1p_CP10K': mean, 'detection_fraction': fraction}))
        evidence[f'sc_{name}_mean_log1p_CP10K'] = np.nanmean(means, axis=0) if means else np.nan
        evidence[f'sc_{name}_mean_detection'] = np.nanmean(fractions, axis=0) if fractions else np.nan
        metadata[name] = {'n_cells': int(mask.sum()), 'n_patients': len(patients)}
    pd.concat(patient_frames).to_csv(OUT / 'single_cell_by_patient.csv', index=False)
    annotation.to_csv(OUT / 'single_cell_annotation.csv', index_label='cell')
    evidence.to_csv(OUT / 'single_cell_gene_evidence.csv', index_label='gene')
    return evidence, {'accession': 'GSE132465', 'n_cells': len(cells), 'n_patients': annotation.Patient.nunique(),
                      'compartments': metadata, 'normalization': 'log1p(UMI / all-gene cell library size * 10000); equal patient averaging',
                      'scope': 'Published cell labels retained; epithelial tumor-tissue location is not a new malignant-cell call or resistance-state measurement.'}


def string_graph(genes):
    mapping = pd.read_csv(SRC / 'STRING_mapping.tsv', sep='\t')
    submitted = list(genes)
    mapping['gene'] = mapping.queryIndex.map(lambda i: submitted[int(i)])
    # Conservative exact HGNC symbol match. Ambiguous aliases cannot silently enter the network.
    accepted = mapping[mapping.preferredName.eq(mapping.gene)].copy()
    lookup = accepted.set_index('stringId').gene.to_dict()
    raw = pd.read_csv(SRC / 'STRING_physical_network.tsv', sep='\t')
    raw['gene_a'] = raw.stringId_A.map(lookup)
    raw['gene_b'] = raw.stringId_B.map(lookup)
    edges = raw.dropna(subset=['gene_a', 'gene_b'])
    edges = edges[edges.gene_a.ne(edges.gene_b) & edges.score.ge(0.7)].copy()
    edges[['gene_a', 'gene_b']] = np.sort(edges[['gene_a', 'gene_b']].to_numpy(), axis=1)
    edges = edges.groupby(['gene_a', 'gene_b'], as_index=False).score.max()
    edges.to_csv(OUT / 'string_edges.csv', index=False)
    mapping.assign(accepted_exact_symbol=mapping.preferredName.eq(mapping.gene)).to_csv(OUT / 'string_mapping_audit.csv', index=False)
    graphs = dict(np.load(OUT / 'graphs.npz'))
    index = {gene: i for i, gene in enumerate(genes)}
    graphs['string'] = edges[['gene_a', 'gene_b']].apply(lambda col: col.map(index)).to_numpy(np.int64).T
    np.savez_compressed(OUT / 'graphs.npz', **graphs)
    return {'version': '12.5', 'network_type': 'physical', 'score_threshold': 0.7,
            'mapped_exact_genes': len(accepted), 'undirected_edges': len(edges),
            'covered_genes': len(set(edges.gene_a) | set(edges.gene_b)),
            'operator': 'Binary thresholded adjacency with self-loops and symmetric degree normalization.'}


def main():
    genes = pd.read_csv(ROOT / 'data/processed/model_genes.csv').gene.tolist()
    annotation = pd.read_csv(OUT / 'genes_annotation.csv').set_index('gene').loc[genes]
    blocks, metadata = [], {}
    for name, prepare in [('depmap', depmap), ('tcga', tcga), ('single_cell', single_cell)]:
        block, details = prepare(genes)
        blocks.append(block)
        metadata[name] = details
        print(json.dumps({name: details}), flush=True)
    metadata['string'] = string_graph(genes)
    for name in ['flux_gene_evidence.csv', 'histone_gene_evidence.csv']:
        blocks.append(pd.read_csv(OUT / name, index_col='gene').reindex(genes))
    metadata.update(json.loads((OUT / 'flux_epigenome_manifest.json').read_text()))
    evidence = pd.concat(blocks, axis=1).reindex(genes)
    evidence.to_csv(OUT / 'integrated_gene_evidence.csv', index_label='gene')
    base = np.load(ROOT / 'data/processed/model_arrays.npz')['gene_features'].copy()
    base[:, -3] = annotation.is_strict_metabolic_enzyme
    base[:, -2] = annotation.is_direct_histone_enzyme
    values = evidence.to_numpy(np.float32)
    features = np.concatenate([base, values, (~np.isfinite(values)).astype(np.float32)], axis=1)
    np.savez_compressed(OUT / 'integrated_features.npz', features=features)
    names = pd.read_csv(ROOT / 'data/processed/model_pathways.csv').pathway_id.tolist() + [
        'is_strict_metabolic_enzyme', 'is_direct_histone_enzyme', 'has_pathway_annotation']
    names += evidence.columns.tolist() + ['missing_' + col for col in evidence.columns]
    pd.DataFrame({'feature_index': range(len(names)), 'feature': names}).to_csv(OUT / 'integrated_feature_names.csv', index=False)
    metadata.update(prepared_utc=datetime.now(timezone.utc).isoformat(), genes=len(genes), features=features.shape[1],
                    normalization='Training-node moments fitted separately per split; missing values standardized to zero with explicit missing indicators.',
                    pair_label_use=False)
    (OUT / 'integrated_input_manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
