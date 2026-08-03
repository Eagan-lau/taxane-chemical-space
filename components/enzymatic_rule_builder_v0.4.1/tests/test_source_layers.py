from enzymatic_rule_builder.source_layers import infer_evidence_layer, infer_evidence_layer_from_record, is_chemical_like_source


def test_bionavi_biochem_is_not_chemical_like():
    assert not is_chemical_like_source("BioNaviNP_BioChem")
    assert infer_evidence_layer("BioNaviNP_BioChem") == "T2_Bio_Extended"


def test_uspto_np_like_is_t3():
    assert is_chemical_like_source("BioNaviNP_USPTO_NPL")
    assert infer_evidence_layer("BioNaviNP_USPTO_NPL") == "T3_Chem_like"


def test_kegg_annotation_only_is_not_promoted_to_t1():
    layer = infer_evidence_layer_from_record({
        "source_database": "KEGG",
        "kegg_ids": "R00001",
        "ec_numbers": "1.1.1.1",
        "reaction_smarts": "",
    })
    assert layer == "T2_Bio_Extended"


def test_retrorules_smarts_with_ec_is_t1_direct_template():
    layer = infer_evidence_layer_from_record({
        "source_database": "RetroRules",
        "source_evidence_text": "retrorules_sqlite_join",
        "reaction_smarts": "[C:1][OH:2]>>[C:1]=[O:2]",
        "template_ec_candidates": "1.1.1.-",
    })
    assert layer == "T1_Bio_Core"


def test_biosynthesis_is_not_chemical_like():
    assert not is_chemical_like_source("taxane biosynthesis pathway")
    assert not is_chemical_like_source("biosynthetic reaction")
