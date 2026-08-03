from pathlib import Path

import pandas as pd

from enzymatic_rule_builder.chem import replay_reaction_smarts_on_pair
from enzymatic_rule_builder.release import export_reaction_smarts_library, validate_smarts_release


def test_replay_validation_accepts_original_pair_direction():
    ok, note = replay_reaction_smarts_on_pair('[C:1][OH:2]>>[C:1]=[O:2]', 'CCO', 'CC=O')
    assert ok is True
    assert note == 'replay_product_match'


def test_smarts_release_has_no_is_reversible_column():
    rules = pd.DataFrame([
        {
            'rule_id': 'RULE_1',
            'reaction_smarts': '[C:1][OH:2]>>[C:1]=[O:2]',
            'template_scope': 'generalized_template',
            'predictive_rule_use': True,
            'template_qc_status': 'ok',
            'strict_core_use': True,
            'expanded_use': True,
            'exploratory_use': True,
            'smarts_direction': 'left_to_right',
            'molecular_direction': 'substrate_to_product',
        }
    ])
    out = export_reaction_smarts_library(rules, tier='core')
    assert len(out) == 1
    assert 'is_reversible' not in out.columns
    assert out.loc[0, 'smarts_direction'] == 'left_to_right'


def test_exact_derived_smarts_without_replay_validation_is_flagged():
    df = pd.DataFrame([
        {
            'rule_id': 'RULE_1',
            'reaction_smarts': '[C:1][OH:2]>>[C:1]=[O:2]',
            'template_scope': 'generalized_template',
            'predictive_rule_use': True,
            'abstracted_from_exact_reaction': True,
            'abstracted_smarts_applies_to_original_pair': False,
            'exact_abstraction_qc_status': 'failed',
        }
    ])
    issues, summary = validate_smarts_release(df)
    assert summary['is_valid_smarts_release'] is False
    assert 'exact_derived_smarts_missing_replay_validation' in set(issues['validation_issue'])


def test_reverse_source_smarts_is_reversed_before_release():
    from enzymatic_rule_builder.templates import build_templates

    df = pd.DataFrame([
        {
            'record_id': 'REC_REV',
            'source_database': 'RetroRules',
            'evidence_layer': 'T1_Bio_Core',
            'source_reaction_id': 'R_REV',
            'reaction_smarts': '[C:1]=[O:2]>>[C:1][OH:2]',
            'direction': 'reverse',
        }
    ])
    templates = build_templates(df)
    assert len(templates) == 1
    assert templates.loc[0, 'reaction_smarts'] == '[C:1][OH:2]>>[C:1]=[O:2]'
    assert templates.loc[0, 'direction_handling'] == 'reversed_from_source'


def test_reversible_source_smarts_is_split_into_two_directional_rules():
    from enzymatic_rule_builder.templates import build_templates, qc_templates
    from enzymatic_rule_builder.rules import build_rule_library

    df = pd.DataFrame([
        {
            'record_id': 'REC_BOTH',
            'source_database': 'RetroRules',
            'evidence_layer': 'T1_Bio_Core',
            'source_reaction_id': 'R_BOTH',
            'reaction_smarts': '[C:1][OH:2]>>[C:1]=[O:2]',
            'direction': 'reversible',
            'reaction_type_source': 'alcohol_oxidoreduction',
            'ec_numbers': '1.1.1.1',
            'template_ec_candidates': '1.1.1.1',
        }
    ])
    templates = qc_templates(build_templates(df))
    assert set(templates['reaction_smarts']) == {'[C:1][OH:2]>>[C:1]=[O:2]', '[C:1]=[O:2]>>[C:1][OH:2]'}
    rules = build_rule_library(templates)
    assert len(rules) == 2
    assert rules['reverse_transform_available'].astype(str).str.lower().eq('true').all()
    assert set(rules['direction_variant']) == {'forward', 'reverse'}
    assert rules['ec_directionality_scope'].eq('direction_specific_rule').all()
    assert rules['ec_directionality_warning'].astype(str).str.contains('must not reuse', case=False).all()
    assert rules['reverse_ec_inheritance_policy'].astype(str).str.contains('do_not_reuse_for_reverse_edge').all()
