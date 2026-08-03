from enzymatic_rule_builder.ec import ec_prefix_match, summarize_ec_evidence


def test_ec_prefix_match_is_hierarchical_not_string_startswith():
    assert ec_prefix_match("1.14", "1.14.14.176")
    assert not ec_prefix_match("1.1", "1.14.14.176")
    assert ec_prefix_match("2.3.1.-", "2.3.1.74")


def test_ec_evidence_ranking_deduplicates_and_prefers_full_supported_ec():
    row = {
        "ec_numbers": "2.3.1.-;2.3.1.74",
        "template_ec_candidates": "2.3.1.74;2.3.1.74",
        "database_ec_candidates": "1.14.14.176",
        "ec_prior_candidates": "2.3.1.-",
        "source_database": "Rhea;RetroRules",
        "evidence_layer": "T1_Bio_Core",
    }
    s = summarize_ec_evidence(row, "T1_Bio_Core")
    assert s["top_ec_number"] == "2.3.1.74"
    assert "2.3.1.74" in s["top3_ec_numbers"]
    assert "2.3.1.-" not in s["top3_ec_numbers"]  # parent partial is redundant once full child exists
    assert s["full_ec_numbers"] == "2.3.1.74;1.14.14.176"
    assert s["partial_ec_numbers"] == "2.3.1.-"
    assert s["broad_ec_class_count"] == 2
    assert s["ec_conflict_level"] == "high"


def test_ec_prior_is_low_confidence_fallback_only():
    row = {
        "ec_numbers": "",
        "template_ec_candidates": "",
        "database_ec_candidates": "",
        "ec_prior_candidates": "1.14.-.-",
        "source_database": "RulePrior",
        "evidence_layer": "Unknown",
    }
    s = summarize_ec_evidence(row, "Unknown")
    assert s["candidate_ec_numbers"] == ""
    assert s["top_ec_number"] == "1.14.-.-"
    assert s["top_ec_assignment_mode"] == "prior_only"
    assert not s["strict_ec_annotation_use"]


def test_high_quality_partial_ec_can_be_top_and_strict():
    row = {
        "ec_numbers": "2.3.1.-",
        "template_ec_candidates": "",
        "database_ec_candidates": "2.3.1.74",
        "ec_prior_candidates": "",
        "source_database": "TaxolKnownPathway_Curated",
        "evidence_layer": "T1_Bio_Core",
        "curated_taxol_anchor": "true",
    }
    s = summarize_ec_evidence(row, "T1_Bio_Core")
    assert s["top_ec_number"] == "2.3.1.-"
    assert float(s["top_ec_confidence"]) >= 0.95
    assert s["top_ec_granularity"] == "subsubclass_level"
    assert s["strict_ec_annotation_use"]
    assert s["ec_conflict_level"] == "none"


def test_ec_confidence_is_not_penalized_by_partial_status():
    full = summarize_ec_evidence({"ec_numbers": "2.3.1.74", "source_database": "Rhea", "evidence_layer": "T1_Bio_Core"}, "T1_Bio_Core")
    partial = summarize_ec_evidence({"ec_numbers": "2.3.1.-", "source_database": "Rhea", "evidence_layer": "T1_Bio_Core"}, "T1_Bio_Core")
    assert full["top_ec_confidence"] == partial["top_ec_confidence"]
    assert full["top_ec_granularity"] == "full4"
    assert partial["top_ec_granularity"] == "subsubclass_level"


def test_reaction_type_ec_class_consistency_blocks_strict_when_incompatible():
    row = {
        "ec_numbers": "1.14.14.176",
        "source_database": "Rhea",
        "evidence_layer": "T1_Bio_Core",
        "reaction_type": "acetylation_or_deacetylation_like_acyl_transfer",
    }
    s = summarize_ec_evidence(row, "T1_Bio_Core")
    assert s["top_ec_number"] == "1.14.14.176"
    assert s["ec_reaction_type_expected_classes"] == "EC2_transferase"
    assert s["ec_reaction_type_observed_classes"] == "EC1_oxidoreductase"
    assert s["ec_reaction_type_consistency"] == "inconsistent"
    assert not s["strict_ec_annotation_use"]


def test_reaction_type_ec_class_consistency_allows_directional_hydrolase_deacetylation():
    row = {
        "ec_numbers": "3.1.1.-",
        "source_database": "Rhea",
        "evidence_layer": "T1_Bio_Core",
        "reaction_type": "deacetylation_or_acetyl_ester_hydrolysis",
    }
    s = summarize_ec_evidence(row, "T1_Bio_Core")
    assert s["ec_reaction_type_expected_classes"] == "EC3_hydrolase"
    assert s["ec_reaction_type_observed_classes"] == "EC3_hydrolase"
    assert s["ec_reaction_type_consistency"] == "consistent"
    assert s["strict_ec_annotation_use"]
