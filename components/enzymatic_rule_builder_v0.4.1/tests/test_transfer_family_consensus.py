import pandas as pd

from enzymatic_rule_builder.chem import reaction_delta, replay_reaction_smarts_on_pair
from enzymatic_rule_builder.transfer_family import build_transfer_family_consensus_templates


def test_transfer_family_consensus_extracts_replayable_o_xylosylation():
    substrate = "CCO"
    product = "CCO[C@@H]1OC[C@@H](O)[C@H](O)[C@H]1O"
    fp, js = reaction_delta(substrate, product)
    main_pairs = pd.DataFrame(
        [
            {
                "record_id": "toy_udp_xyl_1",
                "source_database": "UnitTestBioDB",
                "evidence_layer": "T1_Bio_Core",
                "source_reaction_id": "RXN-XYL",
                "reaction_smiles": substrate + ">>" + product,
                "main_substrate_smiles": substrate,
                "main_product_smiles": product,
                "canonical_reaction_smiles": substrate + ">>" + product,
                "canonical_substrate_smiles": substrate,
                "canonical_product_smiles": product,
                "reaction_delta_fingerprint": fp,
                "reaction_delta_json": js,
                "donor_class": "UDP_sugar_or_nucleotide_sugar",
                "cofactor_or_donor_class": "UDP_sugar_or_nucleotide_sugar",
                "external_participant_roles": "UDP_sugar_or_nucleotide_sugar:left_required_external_participant;UDP_or_nucleotide_diphosphate:right_external_product",
                "transferred_group": "glycosyl",
                "transferred_group_class": "glycosyl",
                "leaving_group_class": "UDP_or_nucleotide_diphosphate",
                "acceptor_atom_class": "alcohol_O_acceptor",
                "main_pair_projection_method": "role_aware_nucleotide_sugar_acceptor_projection",
                "direction": "forward",
                "direction_handling": "kept_forward",
                "direction_variant": "forward",
                "normalized_direction": "substrate_to_product",
                "direction_qc_status": "direction_qc_ok",
            }
        ]
    )

    templates, report, summary = build_transfer_family_consensus_templates(main_pairs)

    assert summary["released_template_rows"] == 1
    assert report.loc[0, "extraction_status"] == "released"
    assert templates.loc[0, "reaction_type_source"] == "O_glycosylation"
    assert templates.loc[0, "reaction_subtype_source"] == "O_pentosylation_or_xylosylation"
    assert templates.loc[0, "consensus_qc_status"] == "passed_transfer_family_replay_validation"
    ok, note = replay_reaction_smarts_on_pair(templates.loc[0, "reaction_smarts"], substrate, product)
    assert ok, note
