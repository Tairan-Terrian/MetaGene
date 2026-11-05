#!/usr/bin/env python3
"""Build the primary metabolic-enzyme × histone-modifying-enzyme scope.

Inputs are the cached human Reactome/HGNC/EpiFactors/UniProt/GO snapshots and
the metabolic/chromatin_adjudication.csv manually reviewed decision tables. This script
only annotates the existing 898-gene model universe; it does not claim to be a
genome-wide catalogue. Run: python scripts/annotate_strict.py
"""
from collections import defaultdict
from datetime import datetime, timezone
from functools import lru_cache
import json
from pathlib import Path
import re

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/evidence'
SOURCE = OUT / 'sources'
ANN = ROOT / 'data/annotations'

DEFINITIONS = {
 'not_in_reviewed_metabolic_universe': 'Outside the reviewed Reactome metabolism subset; no absence-of-function claim.'
}


def joined(values):
    return ' || '.join(dict.fromkeys(str(v) for v in values if str(v) and str(v) != '#'))


def go_ontology(path):
    parents = defaultdict(set)
    names, alternate = {}, {}
    current = None
    version = ''
    with path.open() as stream:
        for raw in stream:
            line = raw.rstrip()
            if line.startswith('data-version: '):
                version = line[14:]
            elif line == '[Term]':
                current = None
            elif line.startswith('id: GO:'):
                current = line[4:]
            elif current and line.startswith('name: '):
                names[current] = line[6:]
            elif current and line.startswith('is_a: '):
                parents[current].add(line[6:16])
            elif current and line.startswith('alt_id: '):
                alternate[line[8:]] = current

    @lru_cache(maxsize=None)
    def is_catalytic(go_id):
        go_id = alternate.get(go_id, go_id)
        return go_id == 'GO:0003824' or any(is_catalytic(parent) for parent in parents[go_id])

    return version, names, is_catalytic


def metabolic_hierarchy():
    roots = json.loads((SOURCE / 'reactome_human_hierarchy.json').read_text())
    root = next(node for node in roots if node['stId'] == 'R-HSA-1430728')
    ids = set()

    def walk(node):
        ids.add(node['stId'])
        for child in node.get('children', []):
            walk(child)

    walk(root)
    return ids


def chromatin_role(functions, gene):
    roles = []
    if 'Histone chaperone' in functions:
        roles.append('histone_chaperone')
    if 'Histone modification read' in functions:
        roles.append('histone_reader')
    if 'Chromatin remodeling' in functions or gene in {'ARID1B', 'SMARCA2', 'SMARCA4'}:
        roles.append('chromatin_remodeling_or_complex_subunit')
    if 'cofactor' in functions:
        roles.append('cofactor')
    if 'DNA modification' in functions:
        roles.append('dna_modification')
    if 'Polycomb' in functions:
        roles.append('polycomb_associated')
    if 'Histone modification' in functions and not roles:
        roles.append('histone_modification_associated')
    return ';'.join(roles)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / 'data/processed/model_genes.csv').fillna('')
    adjudication = pd.read_csv(OUT / 'metabolic_adjudication.csv').fillna('').set_index('gene')
    chromatin_decisions = pd.read_csv(OUT / 'chromatin_adjudication.csv').fillna('').set_index('gene')
    accessions = pd.read_csv(SOURCE / 'hgnc_uniprot_query.csv').fillna('').set_index('gene').uniprot_ids.to_dict()
    uniprot = pd.read_csv(SOURCE / 'uniprot_reviewed_898.tsv', sep='\t', dtype=str).fillna('').set_index('Entry')
    epifactors = pd.read_csv(ANN / 'epifactors_curated.csv', dtype=str).fillna('')
    epi_groups = {gene: group for gene, group in epifactors.groupby('gene')}
    membership = pd.read_csv(ANN / 'pathway_membership.csv', dtype=str)
    metabolic_ids = metabolic_hierarchy()
    metabolic_membership = membership[membership.pathway_id.isin(metabolic_ids)]
    pathway_groups = {gene: group for gene, group in metabolic_membership.groupby('gene')}
    go_version, go_names, is_catalytic = go_ontology(SOURCE / 'go-basic.obo')
    singles = pd.read_csv(ROOT / 'data/processed/singles.csv')
    single_groups = {gene: group for gene, group in singles.groupby('gene')}
    rows = []
    for gene_row in genes.to_dict('records'):
        gene = gene_row['gene']
        row = dict(gene_row)
        row['is_broad_metabolic_broad'] = int(row.pop('is_metabolic'))
        row['is_broad_epigenetic_broad'] = int(row.pop('is_epigenetic'))
        accs = [a for a in accessions.get(gene, '').split('|') if a in uniprot.index]
        entries = uniprot.loc[accs]
        row['reviewed_accessions'] = ';'.join(accs)
        row['reviewed_annotation_status'] = 'resolved_reviewed' if accs else 'unresolved_no_reviewed_record'
        row['uniprot_primary_symbols'] = joined(entries['Gene Names (primary)'])
        for source_col, col in [('Protein names', 'protein_names'), ('EC number', 'ec_numbers'),
                                ('Catalytic activity', 'catalytic_activity'),
                                ('Gene Ontology (molecular function)', 'go_molecular_function'),
                                ('Gene Ontology (biological process)', 'go_biological_process'),
                                ('Function [CC]', 'reviewed_function')]:
            row[col] = joined(entries[source_col])
        row['reviewed_source_urls'] = ';'.join(f'https://www.uniprot.org/uniprotkb/{a}/entry' for a in accs)
        mf_ids = sorted(set(re.findall(r'GO:\d+', row['go_molecular_function'])))
        row['go_mf_ids'] = ';'.join(mf_ids)
        catalytic_ids = [term for term in mf_ids if is_catalytic(term)]
        row['go_catalytic_mf_ids'] = ';'.join(catalytic_ids)
        row['go_catalytic_mf_names'] = ';'.join(go_names.get(term, term) for term in catalytic_ids)
        row['has_go_catalytic_mf_evidence'] = int(bool(catalytic_ids))
        if gene in adjudication.index:
            decision = adjudication.loc[gene]
            for field in ['metabolic_class', 'metabolic_subclass', 'metabolic_reason']:
                row[field] = decision[field]
            row['metabolic_review_status'] = 'manually_reviewed_function_and_reaction'
        else:
            row['metabolic_class'] = 'not_in_reviewed_metabolic_universe'
            row['metabolic_subclass'] = ''
            row['metabolic_reason'] = DEFINITIONS['not_in_reviewed_metabolic_universe']
            row['metabolic_review_status'] = 'outside_broad_metabolic_candidate_universe'
        row['is_strict_metabolic_enzyme'] = int(row['metabolic_class'] == 'strict_core_metabolic_enzyme')
        row['is_extended_metabolic_signaling'] = int(row['metabolic_class'] == 'extended_metabolic_signaling')
        row['is_verified_metabolite_catalyst_including_signaling_lipids'] = int(
            row['is_strict_metabolic_enzyme'] or row['metabolic_subclass'] == 'inositol_or_phosphoinositide_signaling_enzyme')
        pg = pathway_groups.get(gene, membership.iloc[:0])
        row['metabolism_descendant_pathway_ids'] = ';'.join(pg.pathway_id)
        row['metabolism_descendant_pathway_names'] = joined(pg.pathway_name)
        row['metabolism_descendant_pathway_count'] = len(pg)
        ef = epi_groups.get(gene, epifactors.iloc[:0])
        for source_col, col in [('Function', 'epifactors_functions'), ('Modification', 'epifactors_modifications'),
                                ('Target', 'epifactors_targets'), ('Specific_target', 'epifactors_specific_targets'),
                                ('PMID_function', 'epifactors_function_pmids'), ('PMID_target', 'epifactors_target_pmids'),
                                ('Comment', 'epifactors_comments')]:
            row[col] = joined(ef[source_col])
        function = row['epifactors_functions']
        chromatin = bool(re.search('Histone|Chromatin|DNA modification|Polycomb', function))
        rna = 'RNA modification' in function
        row['epifactors_source_url'] = 'https://epifactors.autosome.org/public_data/v2.1/EpiGenes_main.csv' if len(ef) else ''
        row['chromatin_annotation_caveat'] = ''
        if gene == 'NAT10':
            chromatin = False
            epigenetic_class = 'rna_and_uncertain_histone_activity'
            reason = 'EpiFactors contains both histone and RNA entries. Reviewed UniProt establishes RNA acetyltransferase activity and reports histone activity in vitro with uncertain in vivo significance. Excluded from conservative chromatin arm.'
        elif chromatin:
            epigenetic_class = 'strict_chromatin_regulator'
            reason = 'Curated EpiFactors histone/chromatin/DNA-modification/Polycomb function; noncatalytic chaperones, readers and cofactors are eligible as regulators.'
        elif rna:
            epigenetic_class = 'rna_modification_only'
            reason = 'EpiFactors function is RNA modification without a qualifying chromatin function; excluded from strict chromatin partner set.'
        elif len(ef):
            epigenetic_class = 'other_epifactors_not_strict_chromatin'
            reason = 'EpiFactors record has no qualifying histone/chromatin/DNA-modification/Polycomb function; generic TF/scaffold membership alone is insufficient.'
        else:
            epigenetic_class = 'not_in_epifactors_catalogue'
            reason = 'No matched EpiFactors 2.1 main-protein record; absence is not proof of no epigenetic function.'
        caveats = {
            'ARID1B': 'EpiFactors says histone write; reviewed UniProt identifies a SWI/SNF chromatin remodeling subunit. Retained as chromatin regulator, not inferred to be an enzyme.',
            'BRCA2': 'EpiFactors histone write label does not establish direct enzymatic activity; reviewed UniProt primarily describes homologous recombination/DNA repair. Regulator classification follows curated EpiFactors and warrants mechanistic review.',
            'KDM8': 'EpiFactors histone-demethylase label conflicts with current reviewed UniProt histone-tail protease/protein hydroxylase interpretation; retained as histone-associated regulator without endorsing demethylase activity.',
            'PKM': 'Chromatin cofactor annotation concerns moonlighting PKM2 isoform. Gene-level perturbation does not isolate this activity from canonical pyruvate kinase catalysis.',
            'PELP1': 'Histone reader/chromatin cofactor; no direct metabolic enzyme claim. Evidence: EpiFactors reader/write-cofactor and reviewed UniProt chromatin-binding function.',
            'NAT10': reason if gene == 'NAT10' else '',
        }
        row['chromatin_annotation_caveat'] = caveats.get(gene, '')
        row['epigenetic_class'] = epigenetic_class
        row['epigenetic_reason'] = reason
        row['is_chromatin_regulator'] = int(chromatin)
        row['is_rna_modification_factor'] = int(rna)
        row['chromatin_role'] = chromatin_role(function, gene) if chromatin else ''
        if gene in chromatin_decisions.index:
            cd = chromatin_decisions.loc[gene]
            for field in ['is_core_chromatin_factor', 'chromatin_evidence_tier', 'chromatin_mechanism_class',
                          'chromatin_tier_reason', 'chromatin_tier_source', 'chromatin_tier_function_pmids']:
                row[field] = cd[field]
        else:
            row['is_core_chromatin_factor'] = 0
            row['chromatin_evidence_tier'] = 'not_in_epifactors'
            row['chromatin_mechanism_class'] = ''
            row['chromatin_tier_reason'] = 'Outside the manually reviewed EpiFactors candidate universe; no absence-of-function claim.'
            row['chromatin_tier_source'] = ''
            row['chromatin_tier_function_pmids'] = ''
        sg = single_groups.get(gene, singles.iloc[:0])
        row['single_perturbation_datasets'] = ';'.join(sorted(sg.dataset.unique()))
        row['has_ht29_single_perturbation'] = int('ht29_colo1' in set(sg.dataset))
        rows.append(row)
    annotation = pd.DataFrame(rows)
    direct = {
        'histone_acetyltransferase', 'histone_deacetylase', 'histone_demethylase',
        'histone_arginine_methyltransferase', 'histone_ubiquitin_ligase',
        'histone_deubiquitinase', 'histone_acetyltransferase_and_transcription_scaffold',
    }
    annotation['is_direct_histone_enzyme'] = annotation.chromatin_mechanism_class.isin(direct).astype(int)
    annotation['histone_scope_reason'] = annotation.chromatin_mechanism_class.map(
        lambda x: 'Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope.'
        if x in direct else 'Not admitted to the conservative direct histone enzyme set; see preserved mechanism and function evidence.')
    annotation.loc[annotation.gene.eq('KDM8'), 'histone_scope_reason'] = (
        'Mechanism dispute: reviewed histone-tail protease/protein hydroxylase; excluded from the classical histone writer/eraser primary set.')
    annotation.to_csv(OUT / 'genes_annotation.csv', index=False)
    annotation[annotation.is_direct_histone_enzyme.eq(1)].to_csv(OUT / 'direct_histone_enzyme_audit.csv', index=False)
    annotation[annotation.reviewed_annotation_status != 'resolved_reviewed'].to_csv(OUT / 'unresolved_gene_annotations.csv', index=False)
    ht29 = set(singles.loc[singles.dataset.eq('ht29_colo1'), 'gene'])
    sets = {col: set(annotation.loc[annotation[col].eq(1), 'gene']) for col in [
        'is_strict_metabolic_enzyme', 'is_direct_histone_enzyme', 'is_chromatin_regulator',
        'is_core_chromatin_factor', 'is_extended_metabolic_signaling']}
    coverage = {'n_model_genes': len(annotation),
        'model_gene_counts': {k: len(v) for k, v in sets.items()},
        'ht29_single_gene_counts': {k: len(v & ht29) for k, v in sets.items()},
        'direct_histone_genes': sorted(sets['is_direct_histone_enzyme']),
        'unresolved_model_symbols': annotation.loc[annotation.reviewed_annotation_status.ne('resolved_reviewed'), 'gene'].tolist(),
        'scope': 'Existing measured-pair gene universe only; not a genome-wide catalogue.'}
    (OUT / 'annotation_coverage.json').write_text(json.dumps(coverage, indent=2) + '\n')
    definitions = {
        'primary': 'Direct metabolic catalyst x direct histone-modifying catalytic enzyme. Requires individually reviewed mechanism, not a broad pathway or complex label.',
        'epigenetic_extension': 'Metabolic enzyme x other chromatin/epigenetic regulator, including readers, chaperones, remodelers and DNA methyltransferases. Never included in primary Top500.',
        'signaling_extension': 'Extended metabolic signaling regulator x chromatin/epigenetic factor; separate exploratory output.',
        'KDM8': 'Mechanism caveat retained; excluded from primary histone writer/eraser scope.',
        'universe': '38 strictly reviewed metabolic enzymes and the directly annotated histone enzymes among the existing model genes. No missing annotation is interpreted as absent biology.',
        'label_independence': True,
    }
    (OUT / 'annotation_definitions.json').write_text(json.dumps(definitions, indent=2) + '\n')
    manifest = {'prepared_utc': datetime.now(timezone.utc).isoformat(),
        'manual_decisions': ['metabolic_adjudication.csv', 'chromatin_adjudication.csv'],
        'uniprot_downloads': json.loads((SOURCE / 'uniprot_download_manifest.json').read_text()),
        'go_ontology_version': go_version,
        'source_manifest': '../annotations/source_manifest.json',
        'primary_enzymes_audit': 'direct_histone_enzyme_audit.csv',
        'classification_inputs_independent_of_pair_labels': True}
    (OUT / 'annotation_source_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(coverage, indent=2))


if __name__ == '__main__':
    main()
