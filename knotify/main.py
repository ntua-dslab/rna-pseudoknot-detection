import argparse
from datetime import datetime, timedelta
from typing import Tuple

import pandas as pd

from knotify import rna_analysis
from knotify.energy import apply_free_energy_and_stems_criterion


def get_results(
    sequence: str,
    grammar: str = None,
    allow_ug: bool = False,
    max_loop_size: int = rna_analysis.MAX_LOOP_SIZE,
    save_csv: str = None,
    max_stem_allow_smaller: int = 1,
    prune_early: bool = False,
) -> pd.DataFrame:
    """
    Analyze RNA sequence, and predict structure. Return data frame of results
    """
    sequence = sequence.lower()
    pseudoknots = (
        rna_analysis.StringAnalyser(
            input_string=sequence,
            grammar=grammar,
            max_loop_size=max_loop_size,
            allow_ug=allow_ug,
        )
        .get_window_boundaries()
        .generate_trees_in_parallel()
        .get_pseudoknots(
            max_stem_allow_smaller=max_stem_allow_smaller,
            prune_early=prune_early,
        )
    )

    # TODO(akolaitis): consider embedding this into get_pseudoknots()
    # to benefit from multiprocessing library.
    knot_dict_list = [knot.to_dict() for knot in pseudoknots]

    data = pd.DataFrame(knot_dict_list)
    if save_csv is not None:
        data.to_csv(save_csv)

    data = apply_free_energy_and_stems_criterion(
        data, sequence, max_stem_allow_smaller=max_stem_allow_smaller
    )

    return data


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    # store intermediate results
    parser.add_argument("--csv")
    parser.add_argument("--results-csv")

    # pseudoknot arguments
    parser.add_argument("--grammar")
    parser.add_argument("--allow-ug", default=False, action="store_true")
    parser.add_argument("--max-loop-size", default=100, type=int)
    parser.add_argument("--max-stem-allow-smaller", default=2, type=int)
    parser.add_argument("--prune-early", default=False, action="store_true")

    return parser


def main():

    # define the argument parser
    parser = argument_parser()
    parser.add_argument("sequence")

    args = parser.parse_args()

    start = datetime.now()
    results = get_results(
        sequence=args.sequence.lower(),
        grammar=args.grammar,
        save_csv=args.csv,
        allow_ug=args.allow_ug,
        max_loop_size=args.max_loop_size,
        max_stem_allow_smaller=args.max_stem_allow_smaller,
        prune_early=args.prune_early,
    )
    duration = datetime.now() - start

    if args.results_csv:
        results[["dot_bracket", "stems", "energy"]].to_csv(args.results_csv, index=None)

    chosen = results.loc[0]

    print("Sequence: ", args.sequence)
    print("Structure:", chosen.dot_bracket)
    print("Energy:", chosen.energy)
    print("Duration:", duration.total_seconds(), "s")


if __name__ == "__main__":
    main()
