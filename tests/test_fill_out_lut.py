# -*- coding: utf-8 -*-
#
# Tests for the rewritten `cana.utils.fill_out_lut`.
#
# Two things are checked:
#   1. Equivalence: for every input WITHOUT the '2' wildcard, the new (fast)
#      implementation must be byte-for-byte identical to the original one. A
#      verbatim copy of the old implementation is kept below as the reference.
#   2. The '2' fix: '2' used to be accepted but never expanded (leaving spurious
#      non-binary keys and dropped coverage). It must now behave exactly like the
#      other wildcards ('-', '#', 'x').
import random

import pytest

from cana.utils import fill_out_lut


def _reference_fill_out_lut(partial_lut, fill_clashes=False, verbose=False):
    """Verbatim copy of the ORIGINAL fill_out_lut, used as the equivalence oracle.

    Kept intentionally unchanged (including its quirks) so the new implementation
    can be proven equivalent on all non-'2' inputs.
    """
    if len(set([len(x[0]) for x in partial_lut])) != 1:
        raise ValueError(
            "All the input entries of the partial LUT must be of the same length."
        )

    k = len(partial_lut[0][0])
    all_states = dict(partial_lut)

    for entry in partial_lut:
        if not all([x in ["0", "1", "-", "#", "2", "x"] for x in entry[0]]):
            raise ValueError(
                "All the input entries of the partial LUT must be valid binary strings."
            )
        elif any([x in ["-", "#", "2", "x"] for x in entry[0]]):
            missing_data_indices = [
                i for i, x in enumerate(entry[0]) if x in ["-", "#", "x"]
            ]
            table = []
            output_list_permutations = []
            for i in range(2 ** len(missing_data_indices)):
                row = [int(x) for x in bin(i)[2:].zfill(len(missing_data_indices))]
                table.append(row)
                output_list_permutations.append(entry[0])
                for j in range(len(missing_data_indices)):
                    output_list_permutations[i] = (
                        output_list_permutations[i][: missing_data_indices[j]]
                        + str(table[i][j])
                        + output_list_permutations[i][missing_data_indices[j] + 1 :]
                    )
            if entry[0] in all_states:
                del all_states[entry[0]]
            for perm in output_list_permutations:
                if perm in all_states and all_states[perm] != entry[1]:
                    if fill_clashes is False:
                        all_states[perm] = "!"
                    elif fill_clashes is True:
                        all_states[perm] = entry[1]
                else:
                    all_states[perm] = entry[1]

    for i in range(2 ** k):
        state = bin(i)[2:].zfill(k)
        if state not in all_states:
            all_states[state] = "?"

    return sorted(all_states.items(), key=lambda x: x[0])


def _random_partial_lut(rng, k, wildcards):
    """Build a random partial LUT of input length k using the given wildcard chars.

    Produces a mix of concrete and wildcard patterns, and deliberately allows
    repeated / overlapping patterns so clash handling is exercised.
    """
    n_entries = rng.randint(1, 8)
    lut = []
    for _ in range(n_entries):
        pattern = []
        for _ in range(k):
            r = rng.random()
            if r < 0.3 and wildcards:
                pattern.append(rng.choice(wildcards))
            else:
                pattern.append(rng.choice(["0", "1"]))
        value = rng.choice(["0", "1"])
        lut.append(("".join(pattern), value))
    return lut


@pytest.mark.parametrize("k", [2, 3, 4, 5, 6, 7, 8])
@pytest.mark.parametrize("fill_clashes", [False, True])
def test_equivalence_with_reference(k, fill_clashes):
    """New implementation == original on all inputs without the '2' wildcard."""
    rng = random.Random(1234 + k + (1000 if fill_clashes else 0))
    for _ in range(200):
        lut = _random_partial_lut(rng, k, wildcards=["-", "#", "x"])
        got = fill_out_lut(lut, fill_clashes=fill_clashes)
        expected = _reference_fill_out_lut(lut, fill_clashes=fill_clashes)
        assert got == expected, (
            "Mismatch for lut=%s fill_clashes=%s:\n new=%s\n ref=%s"
            % (lut, fill_clashes, got, expected)
        )


def test_documented_example():
    assert fill_out_lut([("00", "0"), ("01", "0"), ("1-", "1"), ("11", "1")]) == [
        ("00", "0"),
        ("01", "0"),
        ("10", "1"),
        ("11", "1"),
    ]


def test_missing_and_clash_markers():
    # '10','11' unspecified -> '?'; a real clash -> '!'
    assert fill_out_lut([("00", "1"), ("01", "1")]) == [
        ("00", "1"),
        ("01", "1"),
        ("10", "?"),
        ("11", "?"),
    ]
    clashed = fill_out_lut([("00--", "0"), ("1--1", "1"), ("11--", "0")])
    assert dict(clashed)["1101"] == "!"
    assert dict(clashed)["1111"] == "!"


@pytest.mark.parametrize("k", [2, 3, 4, 5])
def test_two_wildcard_now_expands_like_dash(k):
    """'2' must now behave identically to '-' (the fix)."""
    rng = random.Random(99 + k)
    for _ in range(100):
        base = _random_partial_lut(rng, k, wildcards=["-"])
        # Mirror the same LUT but swap every '-' for a '2'.
        swapped = [(p.replace("-", "2"), v) for p, v in base]
        assert fill_out_lut(swapped) == fill_out_lut(base)


def test_two_wildcard_no_junk_keys():
    """A pattern whose only wildcard is '2' must expand, not leak a '12' key."""
    result = fill_out_lut([("12", "1")])
    # Every key is a proper length-2 binary string (no '2'/junk keys).
    assert all(set(state) <= {"0", "1"} for state, _ in result)
    assert result == [("00", "?"), ("01", "?"), ("10", "1"), ("11", "1")]
    # Fully-wildcard '2' pattern fills the whole space.
    assert fill_out_lut([("22", "0")]) == [
        ("00", "0"),
        ("01", "0"),
        ("10", "0"),
        ("11", "0"),
    ]


def test_length_and_char_validation():
    with pytest.raises(ValueError):
        fill_out_lut([("00", "0"), ("011", "1")])  # mismatched lengths
    with pytest.raises(ValueError):
        fill_out_lut([("0a", "1")])  # invalid character
