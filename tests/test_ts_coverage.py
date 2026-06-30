# -*- coding: utf-8 -*-
#
# Regression tests for the optimized `computes_ts_coverage`.
#
# These lock the rewritten (loop-inverted) implementation to the exact output
# of the original per-state algorithm: same dict keys, same key order, and the
# same per-state list of covering two-symbol schemata (including order).
#
# The reference implementation below is a verbatim copy of the original
# `computes_ts_coverage` body and relies on `_ts_covers` (which is intentionally
# left unchanged), so it expands schemata the original way. Equality between the
# reference and the optimized function is therefore a byte-for-byte check.
#
# NOTE: these tests require the `schematodes` package (used by
# `find_two_symbols_v2`) to be installed.
import math

import pytest

from cana.boolean_node import BooleanNode
from cana.canalization.boolean_canalization import (
    _ts_covers,
    computes_ts_coverage,
)
from cana.cutils import statenum_to_binstate


def _reference_ts_coverage(k, outputs, two_symbols):
    """Verbatim copy of the original per-state `computes_ts_coverage`."""
    ts_coverage = {}
    for statenum in range(2 ** k):
        binstate = statenum_to_binstate(statenum, base=k)
        ts_coverage[binstate] = covering_twosymbols = []
        output = int(outputs[statenum])
        if output == 2:
            output = [0, 1]
        else:
            output = [int(outputs[statenum])]
        for t in output:
            for implicant, permut_indxs, same_symbols_indxs in two_symbols[t]:
                if _ts_covers(implicant, permut_indxs, binstate):
                    covering_twosymbols.append(
                        (implicant, permut_indxs, same_symbols_indxs)
                    )
    return ts_coverage


def _coverage_pair(outputs):
    """Build a node from an output string and return (new, reference) coverage."""
    k = int(round(math.log(len(outputs), 2)))
    node = BooleanNode(k=k, inputs=list(range(k)), outputs=list(outputs))
    node._check_compute_canalization_variables(two_symbols=True)
    two_symbols = node._two_symbols
    new = computes_ts_coverage(k, node.outputs, two_symbols)
    reference = _reference_ts_coverage(k, node.outputs, two_symbols)
    return new, reference


# A spread of functions across input arities, including permutable groups and
# don't-care symmetries. AP3 (k=7) and Thaliana_AG (k=9) exercise the largest
# expansions; IL7r (k=6) and RULE_110 (k=3) cover smaller cases.
GOLDEN_FUNCTIONS = {
    "AND": "0001",
    "OR": "0111",
    "XOR": "0110",
    "RULE_110": "01110110",
    "Lymphoid_IL7r": (
        "0000000000000000010101010100010000001111000011000101111101001100"
    ),
    "AP3": (
        "0000000000000000000000000000000000000000111111110000000011111111"
        "0000000100000001000100010001000100000001111111110001000111111111"
    ),
}


@pytest.mark.parametrize("name", sorted(GOLDEN_FUNCTIONS))
def test_ts_coverage_matches_reference(name):
    outputs = GOLDEN_FUNCTIONS[name]
    new, reference = _coverage_pair(outputs)

    # Same keys, in the same (ascending state-number) order.
    assert list(new.keys()) == list(reference.keys()), (
        "Key order differs for %s" % name
    )
    # Same per-state covering lists, including order within each state.
    for binstate in reference:
        assert new[binstate] == reference[binstate], (
            "Coverage for state %s differs for %s:\n  new=%s\n  ref=%s"
            % (binstate, name, new[binstate], reference[binstate])
        )
    # Full dict equality as a final catch-all.
    assert new == reference
