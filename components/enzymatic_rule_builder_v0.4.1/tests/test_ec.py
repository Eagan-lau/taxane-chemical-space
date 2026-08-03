from enzymatic_rule_builder.ec import broad_ec_classes, normalize_ec_numbers


def test_partial_ec_numbers_are_preserved():
    assert normalize_ec_numbers("2.3.1.-") == "2.3.1.-"
    assert normalize_ec_numbers("EC 1.14.-.-; 3.1.1.-") == "1.14.-.-;3.1.1.-"


def test_broad_ec_classes_from_partial():
    assert broad_ec_classes("1.14.-.-;2.3.1.-") == "EC1_oxidoreductase;EC2_transferase"
