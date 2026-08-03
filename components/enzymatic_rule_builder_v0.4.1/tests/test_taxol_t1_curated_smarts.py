import pandas as pd

from enzymatic_rule_builder.rules import build_rule_library


def test_taxol_t1_curated_smarts_can_enter_core_with_partial_ec():
    templates = pd.DataFrame([
        {
            'template_hash': 'h1',
            'template_scope': 'generalized_template',
            'predictive_rule_use': 'true',
            'anchor_edge_use': 'true',
            'template_qc_status': 'ok',
            'template_qc_note': '',
            'reaction_smarts': '[C:1][OH:2]>>[C:1]=[O:2]',
            'canonical_reaction_smiles': 'CCO>>CC=O',
            'main_substrate_smiles': 'CCO',
            'main_product_smiles': 'CC=O',
            'reaction_delta_fingerprint': 'H-2|M-2.016',
            'reaction_type_source': 'alcohol_oxidation',
            'reaction_subtype_source': '',
            'direction': 'forward',
            'source_database': 'TaxolKnownPathway_Curated',
            'evidence_layer': 'T1_Bio_Core',
            'source_reaction_id': 'TAXOL_TEST',
            'record_id': 'REC1',
            'ec_numbers': '1.1.-.-',
            'template_ec_candidates': '',
            'database_ec_candidates': '',
            'ec_prior_candidates': '',
            'rhea_ids': '',
            'kegg_ids': '',
            'metanetx_ids': '',
            'cofactor_or_donor_class': '',
            'enzyme_name': 'test enzyme',
            'protein_ids': '',
            'source_evidence_text': 'curated exact reaction generalized and replay validated',
            'source_file': 'taxol_pathway.csv',
            'curated_taxol_anchor': 'true',
            'curated_pathway_name': 'TaxolKnownPathway',
            'curated_pathway_step_id': 'taxol_pathway_1',
            'abstracted_from_exact_reaction': 'true',
            'derived_from_exact_anchor': 'true',
            'rxnmapper_confidence': '0.98',
            'rdchiral_extraction_status': 'ok',
            'abstracted_smarts_applies_to_original_pair': 'true',
            'exact_abstraction_qc_status': 'pass',
            'benchmark_exclusion_flag': 'exclude_from_external_recall',
        }
    ])
    rules = build_rule_library(templates)
    assert len(rules) == 1
    row = rules.iloc[0]
    assert row['evidence_layer_best'] == 'T1_Bio_Core'
    assert row['curated_taxol_anchor'] == 'true'
    assert row['strict_ec_annotation_use'] in {True, 'true'}
    assert row['strict_core_use'] in {True, 'true'}
    assert row['benchmark_exclusion_flag'] == 'exclude_from_external_recall'
    assert 'is_reversible' not in rules.columns
