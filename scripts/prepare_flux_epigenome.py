"""Map published HCT116/OxR model flux bounds and histone peaks to model genes.

These are independent covariates, not double-perturbation labels. Flux ranges
are author model estimates, and peak summaries are descriptive (no replicate
differential-binding test). GPR membership is not a gene-specific flux estimate.
"""
import importlib
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/evidence'
SRC = OUT / 'requested_sources'


class RestrictedFrameReader(pickle.Unpickler):
    """Only the classes used by these fixed, hash-recorded pandas snapshots."""
    allowed = {('pandas.core.frame', 'DataFrame'),
               ('pandas.core.internals.managers', 'BlockManager'),
               ('pandas.core.indexes.base', '_new_Index'),
               ('pandas.core.indexes.base', 'Index'),
               ('numpy.core.multiarray', '_reconstruct'),
               ('numpy', 'ndarray'), ('numpy', 'dtype'), ('builtins', 'slice')}

    def find_class(self, module, name):
        if (module, name) not in self.allowed:
            raise ValueError(f'Unexpected serialized class: {module}.{name}')
        return getattr(importlib.import_module(module), name)


def flux(genes):
    model = loadmat(SRC / 'HCT116_metabolic_model.mat', simplify_cells=True)
    model = next(value for key, value in model.items() if not key.startswith('_'))
    hgnc = pd.read_csv(ROOT / 'data/annotations/sources/hgnc_complete_set.txt', sep='\t', low_memory=False)
    lookup = hgnc.dropna(subset=['ensembl_gene_id']).drop_duplicates('ensembl_gene_id').set_index('ensembl_gene_id').symbol
    reaction_genes = model['rxnGeneMat'].tocoo()
    members = pd.DataFrame({'reaction': np.asarray(model['rxns'])[reaction_genes.row],
                            'ensembl_gene_id': np.asarray(model['genes'])[reaction_genes.col]})
    members['gene'] = members.ensembl_gene_id.map(lookup)
    members = members[members.gene.isin(genes)].drop_duplicates(['reaction', 'gene'])
    frames, growth = [], {}
    for condition in ['Control', 'OxR']:
        with (SRC / f'HCT116_flux_{condition}_clean.p').open('rb') as stream:
            frame = RestrictedFrameReader(stream).load().astype(float)
        if not np.isfinite(frame.to_numpy()).all() or not (frame.minimum <= frame.maximum + 1e-6).all():
            raise ValueError('Invalid author flux bounds')
        growth[condition] = float(frame.loc['HCC_biomass', 'minimum'])
        if growth[condition] <= 0:
            raise ValueError('Non-positive growth normalization denominator')
        frame = frame.rename(columns={c: f'{condition}_{c}_raw' for c in frame.columns})
        for bound in ['minimum', 'maximum']:
            frame[f'{condition}_{bound}_growth_normalized'] = frame[f'{condition}_{bound}_raw'] / growth[condition]
        frame[f'{condition}_log1p_abs_bound'] = np.log1p(frame[[f'{condition}_{bound}_growth_normalized' for bound in ['minimum', 'maximum']]].abs().max(axis=1))
        frames.append(frame)
    reaction = pd.concat(frames, axis=1, join='inner')
    reaction.index.name = 'reaction'
    reaction['OxR_minus_Control_log1p_abs_bound'] = reaction.OxR_log1p_abs_bound - reaction.Control_log1p_abs_bound
    reaction.to_csv(OUT / 'hct116_flux_reaction_bounds.csv')
    mapped = members.merge(reaction.reset_index(), on='reaction')
    mapped.to_csv(OUT / 'hct116_flux_gene_reaction_map.csv', index=False)
    grouped = mapped.groupby('gene')
    evidence = pd.DataFrame(index=genes)
    evidence['flux_HCT116_Control_mean_log1p_abs_bound'] = grouped.Control_log1p_abs_bound.mean()
    evidence['flux_HCT116_OxR_mean_log1p_abs_bound'] = grouped.OxR_log1p_abs_bound.mean()
    evidence['flux_HCT116_OxR_minus_Control_mean_log1p_abs_bound'] = grouped.OxR_minus_Control_log1p_abs_bound.mean()
    evidence['flux_HCT116_mapped_reaction_count'] = grouped.reaction.nunique()
    evidence.to_csv(OUT / 'flux_gene_evidence.csv', index_label='gene')
    return evidence, {'doi': '10.3390/cancers13164130', 'author_repository': 'https://github.com/HAHerrmann/Hct116_DrugRes',
        'commit': '3ba911ee3b931e0909add599beba48bdca9f69a9', 'author_bound_file_variant': 'pTVA_*_clean.p',
        'reaction_count': len(reaction), 'mapped_reactions': mapped.reaction.nunique(), 'covered_genes': mapped.gene.nunique(),
        'biomass_minimum_denominators': growth,
        'normalization': 'Divide each state bound by its HCC_biomass minimum, following author normalize_growth; log1p(max(abs(lower),abs(upper))), then equal mean over linked reactions.',
        'mapping': 'HCT_GEM.mat rxnGeneMat association -> HGNC Ensembl symbols; only exact shared reaction identifiers. GPR AND/OR semantics are not used to simulate knockout.',
        'scope': 'Author thermodynamic model feasible bounds constrained by metabolomics, not measured intracellular fluxes, statistical intervals, or gene-knockout effects. Clean snapshot reanalysis; not a reproduction of every published figure.'}


def epigenome(genes):
    ref = pd.read_csv(SRC / 'hg19_refGene.txt.gz', sep='\t', header=None, usecols=[1, 2, 3, 4, 5, 12],
                      names=['transcript', 'chr', 'strand', 'tx_start', 'tx_end', 'gene'])
    # Exact current symbols only; all annotated transcript TSS windows are unioned.
    ref = ref[ref.gene.isin(genes) & ref.chr.str.match(r'^chr([0-9]+|X|Y)$')].copy()
    ref['tss'] = np.where(ref.strand.eq('+'), ref.tx_start, ref.tx_end - 1)
    ref = ref.drop_duplicates(['gene', 'chr', 'tss'])
    ref['promoter_start'] = np.maximum(0, ref.tss - 2000)
    ref['promoter_end'] = ref.tss + 2001
    ref.to_csv(OUT / 'histone_promoter_mapping.csv', index=False)
    evidence = pd.DataFrame(index=genes)
    evidence['histone_hg19_promoter_mapped'] = np.asarray([gene in set(ref.gene) for gene in genes], dtype=int)
    samples, links = [], []
    for gsm, mark, state in [(4115600, 'H3K4me3', 'HCT116'), (4115601, 'H3K4me3', 'HCT116OxR'),
                             (4115602, 'H3K27ac', 'HCT116'), (4115603, 'H3K27ac', 'HCT116OxR'),
                             (4115604, 'H3K27me3', 'HCT116'), (4115605, 'H3K27me3', 'HCT116OxR')]:
        peaks = pd.read_csv(SRC / f'GSM{gsm}_{mark}_{state}_peaks.txt.gz', sep='\t')
        # MACS XLS has 1-based inclusive endpoints; refGene is 0-based half-open.
        peaks['start0'] = peaks.start - 1
        by_chr = {chrom: part for chrom, part in peaks.groupby('chr')}
        counts, enrichment = {}, {}
        for gene, promoters in ref.groupby('gene'):
            overlap_ids = set()
            for promoter in promoters.itertuples():
                chrom_peaks = by_chr.get(promoter.chr)
                if chrom_peaks is None:
                    continue
                mask = (chrom_peaks.start0 < promoter.promoter_end) & (chrom_peaks.end > promoter.promoter_start)
                overlap_ids.update(chrom_peaks.index[mask])
            matched = peaks.loc[sorted(overlap_ids)]
            counts[gene] = len(matched)
            enrichment[gene] = float(matched.fold_enrichment.max()) if len(matched) else 0.0
            for peak in matched.itertuples():
                links.append({'gene': gene, 'sample': f'GSM{gsm}', 'state': state, 'mark': mark,
                              'peak': peak.name, 'chr': peak.chr, 'start0': peak.start0, 'end': peak.end,
                              'fold_enrichment': peak.fold_enrichment})
        evidence[f'ChIP_{mark}_{state}_log1p_peak_count'] = np.log1p(pd.Series(counts))
        evidence[f'ChIP_{mark}_{state}_log1p_max_fold_enrichment'] = np.log1p(pd.Series(enrichment))
        samples.append({'sample': f'GSM{gsm}', 'mark': mark, 'state': state, 'peaks': len(peaks),
                        'genes_with_promoter_peak': sum(count > 0 for count in counts.values()),
                        'url': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM{gsm}'})
    for mark in ['H3K4me3', 'H3K27ac', 'H3K27me3']:
        evidence[f'ChIP_{mark}_OxR_minus_parental_log1p_max_fold_enrichment'] = evidence[f'ChIP_{mark}_HCT116OxR_log1p_max_fold_enrichment'] - evidence[f'ChIP_{mark}_HCT116_log1p_max_fold_enrichment']
    evidence.to_csv(OUT / 'histone_gene_evidence.csv', index_label='gene')
    pd.DataFrame(links).to_csv(OUT / 'histone_gene_peak_evidence.csv', index=False)
    pd.DataFrame(samples).to_csv(OUT / 'histone_sample_manifest.csv', index=False)
    return evidence, {'accession': 'GSE138654', 'super_series': 'GSE138657', 'genome': 'hg19',
        'promoter_mapped_genes': ref.gene.nunique(), 'samples': samples,
        'summary': 'Union of +/-2000bp windows around all exact-symbol RefGene TSS; count unique overlapping author MACS2 peaks and maximum fold enrichment, log1p transformed.',
        'scope': 'Six deposited processed peak tracks, one per mark/state. No replicate-based differential-binding test; missing called peak is not proof of absent modification. Differences are descriptive, with no cross-sample spike-in calibration. These HCT116/OxR lines are not assumed identical to collaborators or the separate flux study.'}


def main():
    genes = pd.read_csv(ROOT / 'data/processed/model_genes.csv').gene.tolist()
    metadata = {}
    for name, prepare in [('metabolic_flux', flux), ('epigenetic_map', epigenome)]:
        _, metadata[name] = prepare(genes)
    (OUT / 'flux_epigenome_manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
