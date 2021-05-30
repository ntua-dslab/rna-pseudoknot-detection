import os

import pytest

from knotify.main import get_results

PSEUDOKNOT = os.getenv("PSEUDOKNOT_SO", "./libpseudoknot.so")


@pytest.mark.parametrize(
    "name, sequence, result, test_params",
    [
        (
            "baseline",
            "GGGAAAUGGACUGAGCGGCGCCGACCGCCAAACAACCGGCA",
            "..............((((.[[[[.))))........]]]].",
            {},
        ),
        (
            "baseline, ug",
            "ACGUGAAGGCUACGAUAGUGCCAG",
            ".((((..[[[)))).....]]]..",
            {"allow_ug": True},
        ),
        (
            "baseline",
            "GGGAAACGGAGUGCGCGGCACCGUCCGCGGAACAAACGGAGAAGGCAGCU",
            ".............(((((..[[[[)))))......]]]]...........",
            {},
        ),
        (
            "baseline",
            "AUCCUUUUCAGUUGGGCCUUCUGGUGAUGUUUCUGGCCACCCAGGAGGUCCUGAGGAAGAGGUGGACGGCCAGAUUGACU",
            ".............(((((((((((.......[[[[[[[..)))))))))))................]]]]]]]......",
            {},
        ),
        (
            "pick smaller dd for same energy, stems",
            "AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC",
            "..............((((((.....[[[))))))....]]]......",
            {},
        ),
    ],
)
def test_end_to_end(name, sequence, result, test_params):
    config = {
        "grammar": PSEUDOKNOT,
        "prune_early": True,
    }
    config.update(test_params)

    results = get_results(sequence, **config)
    assert results.loc[0].dot_bracket == result
