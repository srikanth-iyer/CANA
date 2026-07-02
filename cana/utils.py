import copy
import math
import operator as op
from functools import reduce
from itertools import takewhile

import networkx as nx
import numpy as np

from cana.cutils import binstate_to_statenum, flip_binstate_bit, statenum_to_binstate


def flip_bitset_in_strstates(strstates, idxs):
    """Flips the binary value for a set of bits in a binary state.

    Args:
        binstate (string) : The binary state to flip.
        idxs (int) : The indexes of the bits to flip.

    Returns:
        (list) : The flipped states

    Example:
        >>> flip_bit_in_strstates('000',[0, 2])
        ['100','001']
    """
    return [flip_binstate_bit(strstates, idx) for idx in idxs]


def flip_binstate_bit_set(binstate, idxs):
    """Flips the binary value for a set of bits in a binary state.

    Args:
        binstate (string) : The binary state to flip.
        idxs (int) : The indexes of the bits to flip.

    Returns:
        (list) : The flipped states
    """
    flipset = []
    if len(idxs) != 0:
        fb = idxs.pop()
        flipset.extend(flip_binstate_bit_set(binstate, copy.copy(idxs)))
        flipset.extend(
            flip_binstate_bit_set(flip_binstate_bit(binstate, fb), copy.copy(idxs))
        )
    else:
        flipset.append(binstate)
    return flipset


def negate_LUT_input(outputs, idx):
    """For a LUT defined by the output list, it negates the input.

    Args:
        outputs (list) : The output list defining the LUT.
        idxs (int) : The indexes of the input to negate.

    Returns:
        (list) : The new output with input idx negated
    """
    k = int(np.log2(len(outputs)))
    return [
        outputs[
            binstate_to_statenum(
                flip_binstate_bit(statenum_to_binstate(istate, k), idx)
            )
        ]
        for istate in range(len(outputs))
    ]


def print_logic_table(outputs):
    """Print Logic Table

    Args:
        outputs (list) : The transition outputs of the function.

    Returns:
        print : a print-out of the logic table.

    Example:
        >>> print_logic_table([0,0,1,1])
        00 : 0
        01 : 0
        10 : 1
        11 : 1

    """
    k = int(math.log(len(outputs)) / math.log(2))
    for statenum in range(2**k):
        print(statenum_to_binstate(statenum, base=k) + " : " + str(outputs[statenum]))


def entropy(prob_vector, logbase=2.0):
    """Calculates the entropy given a probability vector

    Todo:
        This should be calculated using ``scipy.entropy``
    """
    prob_vector = np.array(prob_vector)
    pos_prob_vector = prob_vector[prob_vector > 0]
    return -np.sum(pos_prob_vector * np.log(pos_prob_vector) / np.log(logbase))


def ncr(n, r):
    """Return the combination number.
    The combination of selecting `r` items from `n` iterms, order doesn't matter.

    Args:
        n (int): number of elements in collection
        r (int): length of combination

    Returns:
        int
    """
    r = min(r, n - r)
    if r == 0:
        return 1
    numer = reduce(op.mul, range(n, n - r, -1))
    denom = reduce(op.mul, range(1, r + 1))
    return numer // denom


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """Python 2 doesn't have math.isclose()
    Here is an equivalent function
    Use this to tell whether two float numbers are close enough
    considering using == to compare floats is dangerous!
    2.0*3.3 != 3.0*2.2 in python!

    Args:
        a (float) : the first float number
        b (float) : the second float number
        rel_tol (float) : the relative difference threshold between a and b
        abs_tol (float) : absolute difference threshold. not recommended for float

    Returns:
        bool
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def output_transitions(eval_line, input_list):
    """Returns an output list from combinatorically trying all input values

    Args:
        eval_line (string) : logic or arithmetic line to evaluate
        input_list (list) : list of input variables

    Returns:
        list of all possible output transitions (list)

    Example:
        RAS*=(GRB2 or PLCG1) and not GAP

        .. code-block:: python

            >>> eval_line = "(GRB2 or PLCG1) and not GAP"
            >>> input_list = ['GRB2', 'PLCG1', 'GAP']
            >>> output_transitions(eval_line, input_list)
            000
            001
            010
            011
            100
            101
            110
            111

        A variable is dynamically created for each member of the input list
        and assigned the corresponding value from each trail string.
        The original eval_line is then evaluated with each assignment
        which results in the output list [0, 0, 1, 0, 1, 0, 1, 0]
    """
    total = 2 ** len(input_list)  # Total combinations to try
    output_list = []
    for i in range(total):
        trial_string = statenum_to_binstate(i, len(input_list))
        # Evaluate trial_string by assigning value to each input variable
        for j, input in enumerate(input_list):
            exec(input + "=" + trial_string[j])
        output_list.append(int(eval(eval_line)))

    return output_list


def seq_upto(seq, obj):
    """
    TODO: description
    """
    return takewhile(lambda el: el != obj, iter(seq))


def mindist_from_source(G, source):
    """
    TODO: description
    """
    dag = nx.bfs_tree(G, source)
    dist = {}  # stores [node, distance] pair
    for node in nx.topological_sort(dag):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0] + 1, v) for v in dag.pred[node]]
        if pairs:
            dist[node] = min(pairs)
        else:
            dist[node] = (0, node)

    return dist


def pathlength(p, weights, rule="sum"):
    """Calculate the length of path p, with weighted edges, given the length rule of:

    Ars:
        weights:

        rule (str):
            'sum' - sum of edge weights along path
            'prod' - product of edge weights along path
            'min' - minimum of edge weights along path (weakest-link)
            'max' - maximum of edge weights along path

    TODO: update description
    """
    if rule == "sum":
        return sum(weights[(p[ie], p[ie + 1])] for ie in range(len(p) - 1))
    elif rule == "prod":
        return np.prod([weights[(p[ie], p[ie + 1])] for ie in range(len(p) - 1)])
    elif rule == "min":
        return min(weights[(p[ie], p[ie + 1])] for ie in range(len(p) - 1))
    elif rule == "max":
        return max(weights[(p[ie], p[ie + 1])] for ie in range(len(p) - 1))


def function_monotone(outputs, method="exact", nsamples=100, random_seed=None):
    """
    Determine if a given LUT is monotone.

    Here we test every pair of inputs that are Hamming distance 1. (see Goldreich et al 2000)

    Args:
        outputs (list) : The transition outputs of the function.

        method (str) :
            'exact' - test all pairs of inputs
            TODO: 'random' - sample pairs of inputs

        nsamples (int) : when method=='random', specifies the number of samples.

    Returns:
        (Bool) : True if monotone.

    Example:
        >>> is_monotone(outputs=[0,0,0,1])
    """
    # np.random.seed(random_seed)

    k = int(np.log2(len(outputs)))

    # for all input configurations
    for input_confignum in range(2**k):
        input_configbin = statenum_to_binstate(input_confignum, k)

        # for all input states that are 0
        for idx, state in enumerate(input_configbin):
            if state == "0":
                # we flip the 0 and check for monotone along the edge
                next_confignum = binstate_to_statenum(
                    flip_binstate_bit(input_configbin, idx)
                )

                # if the monotone property fails
                if outputs[input_confignum] > outputs[next_confignum]:
                    return False
    return True


def input_monotone(outputs, input_idx, activation=1):
    """
    Determine if a given input is activating or inhibiting in a given function.

    Args:
        outputs (list) : The transition outputs of the function.

        input_idx (int) : The input to test.

        activation (1 or -1) : Whether to test for activation or inhibition.

    Returns:
        (Bool) : True if monotone with respect to activation or inhibition.

    Example:
        >>> input_monotone([0,1,0,0], 0, activation=1) == False
        >>> input_monotone([0,1,0,0], 0, activation=-1) == True
    """

    k = int(np.log2(len(outputs)))

    if k == 1:
        return True
    else:
        monotone_configs = []
        # for all input configurations
        for input_confignum in range(2 ** (k - 1)):
            other_input_configbin = statenum_to_binstate(input_confignum, k - 1)

            input_confignum_0 = binstate_to_statenum(
                other_input_configbin[:input_idx]
                + "0"
                + other_input_configbin[input_idx:]
            )

            input_confignum_1 = binstate_to_statenum(
                other_input_configbin[:input_idx]
                + "1"
                + other_input_configbin[input_idx:]
            )

            if activation == 1:
                monotone_configs.append(
                    outputs[input_confignum_0] <= outputs[input_confignum_1]
                )
            elif activation == -1:
                monotone_configs.append(
                    outputs[input_confignum_0] >= outputs[input_confignum_1]
                )

        return all(monotone_configs)


# Sentinel marking a state whose output has not been provided (directly or via a
# wildcard expansion). It is a unique object -- deliberately NOT the string '?' --
# so that a partial-LUT entry that *explicitly* carries the value '?' stays
# distinguishable from a genuinely-unspecified state during clash detection.
_LUT_UNSET = object()

# The characters that stand for a "don't care" position in a partial-LUT input.
# One canonical set is used in every check below. Previously three separate
# hard-coded lists disagreed about '2': it was accepted as valid and routed into
# the wildcard branch but never actually expanded, which left spurious non-binary
# keys (e.g. "12") in the output and dropped the coverage it should have added.
# Using a single set makes '2' behave like the other wildcards.
_LUT_WILDCARDS = frozenset("-#x2")
_LUT_VALID_CHARS = frozenset("01-#x2")


def fill_out_lut(partial_lut: list, fill_clashes: bool = False, verbose: bool = False) -> list:
    """Expand a partial look-up table (LUT) into a complete one.

    A *partial* LUT specifies the output for only some input states, and its input
    patterns may contain wildcard characters (``-``, ``#``, ``x`` or ``2``) that each
    stand for "either 0 or 1". This function expands every wildcard pattern into the
    concrete input states it covers and returns the **full** LUT for all ``2**k``
    input states (``k`` = length of the input patterns), sorted in ascending binary
    order.

    States that are never specified are filled with ``'?'``. When two entries assign
    conflicting outputs to the same state it is marked ``'!'`` (a clash), unless
    ``fill_clashes=True``, in which case the later conflicting entry's value is kept.

    Args:
        partial_lut (list): list of ``(input_pattern, output)`` pairs. Every
            ``input_pattern`` is a string of the same length using the characters
            ``0``, ``1`` or a wildcard (``-``, ``#``, ``x``, ``2``). ``output``
            values are stored verbatim (commonly the strings ``'0'`` / ``'1'``).
        fill_clashes (bool): if ``False`` (default) a conflicting state becomes
            ``'!'``; if ``True`` the later conflicting entry's value overwrites the
            earlier one and no ``'!'`` is emitted.
        verbose (bool): if ``True`` print clash notices and whether the resulting LUT
            is complete.

    Returns:
        list: ``(binstate, output)`` pairs for all ``2**k`` states, sorted by
        ``binstate`` (ascending binary order); missing states are ``'?'`` and
        clashes are ``'!'``.

    Raises:
        ValueError: if the input patterns are not all the same length, or contain a
            character other than ``0``, ``1`` or a wildcard.

    Example:
        >>> fill_out_lut([('00', '0'), ('01', '0'), ('1-', '1'), ('11', '1')])
        [('00', '0'), ('01', '0'), ('10', '1'), ('11', '1')]

    # TODO: [SRI] generate LUT from two symbol schemata, with a specified ratio of wildcard symbols
    """
    # Every input pattern must share one length k, which defines the 2**k state
    # space. An empty partial_lut produces an empty length-set, correctly failing here.
    lengths = {len(pattern) for pattern, _ in partial_lut}
    if len(lengths) != 1:
        raise ValueError(
            "All the input entries of the partial LUT must be of the same length."
        )
    k = next(iter(lengths))
    n_states = 1 << k  # == 2**k

    # Validate every character up front so an invalid pattern fails fast.
    for pattern, _ in partial_lut:
        if not all(c in _LUT_VALID_CHARS for c in pattern):
            raise ValueError(
                "All the input entries of the partial LUT must be valid binary strings."
            )

    # Dense table indexed by integer state number; _LUT_UNSET means "not yet set".
    # Indexing by int(binstate, 2) lets us (a) expand wildcards with cheap integer
    # bit math and (b) emit the result already in sorted order -- replacing the old
    # per-bit string slicing and the final O(2**k log 2**k) sort.
    table = [_LUT_UNSET] * n_states

    # Pass 1 -- concrete entries (no wildcards). Placing these first mirrors the old
    # dict(partial_lut) behaviour: a fully-specified state is present before any
    # wildcard expansion runs, and a later duplicate simply overwrites the earlier.
    for pattern, value in partial_lut:
        if not any(c in _LUT_WILDCARDS for c in pattern):
            table[int(pattern, 2)] = value

    # Pass 2 -- wildcard entries, in list order, so clash resolution matches the
    # original (a later conflicting entry wins as '!' -- or as its own value when
    # fill_clashes=True).
    for pattern, value in partial_lut:
        wildcard_positions = [i for i, c in enumerate(pattern) if c in _LUT_WILDCARDS]
        if not wildcard_positions:
            continue  # concrete entry, already handled in pass 1

        # `base` = the state number with every wildcard bit set to 0. Position i of
        # the MSB-first pattern has integer weight 1 << (k - 1 - i).
        base = 0
        for i, c in enumerate(pattern):
            if c == "1":
                base |= 1 << (k - 1 - i)
        # Integer weight contributed by each wildcard position.
        weights = [1 << (k - 1 - i) for i in wildcard_positions]
        m = len(wildcard_positions)

        # Each of the 2**m subsets of the wildcard bits yields one covered state.
        for combo in range(1 << m):
            idx = base
            for j in range(m):
                if combo & (1 << j):
                    idx |= weights[j]

            current = table[idx]
            if current is not _LUT_UNSET and current != value:
                # State already holds a different specified value -> conflict.
                if verbose:
                    print("Clashing output values for state:", format(idx, "0{}b".format(k)))
                table[idx] = value if fill_clashes else "!"
            else:
                table[idx] = value

    # Build the full LUT. Ascending idx == ascending binstate, so no sort is needed;
    # any state still UNSET is reported as the missing marker '?'. bin(idx)[2:].zfill
    # is the fastest way to render the fixed-width binstate for every state.
    filled_lut = [
        (bin(idx)[2:].zfill(k), "?" if table[idx] is _LUT_UNSET else table[idx])
        for idx in range(n_states)
    ]

    if verbose:
        # A LUT is incomplete if any state ended up as the missing marker '?'.
        if any(value == "?" for _, value in filled_lut):
            print("The LUT is incomplete. Missing values are represented by '?'")
        else:
            print("The LUT is complete.")

    return filled_lut
