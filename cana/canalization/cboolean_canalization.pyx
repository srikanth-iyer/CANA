# -*- coding: utf-8 -*-
"""
(Cythonized) Boolean Canalization
=====================

Functions to compute the Quine-McCluskey algorithm in cython for increaed computation speed.

"""
#   Copyright (C) 2021 by
#   Rion Brattig Correia <rionbr@gmail.com>
#   Alex Gates <ajgates42@gmail.com>
#
#   All rights reserved.
#   MIT license.
import cython

from cana.cutils import *

WILDCARD_SYMBOL = '#'
SYMMETRIC_WILDCARD_SYMBOL = '*'


#
# Quine-McCluskey Functions
#
def make_density_groups(input_binstates):
    """
    """

    density_groups = dict()
    for binstate in input_binstates:
        density = binary_density(binstate)
        if density not in density_groups:
            density_groups[density] = set()
        density_groups[density].add(binstate)

    return density_groups


def find_wildcards(binstate1, binstate2):
    """
    Compare two binary states and replace any differing bits by a wildcard.
    Args:
        binstate1, binstate2 : the two binary states to be compared

    Return:
        c (list, bool) : a list of comparisons

    """
    # assert len(s1) == len(s2) , "The two binstates must have the same length"
    return "".join([b0 if (b0 == b1) else WILDCARD_SYMBOL for b0, b1 in zip(binstate1, binstate2)])


def binary_density(binstate):
    """
    Find the density (number of 1s) for a term with possible wildcards.
    """
    return binstate.count('1')


def replace_wildcard(binstate, idx):
    """
    Return the binary state with a wildcard at the idx position.
    """
    return binstate[:idx] + WILDCARD_SYMBOL + binstate[idx + 1:]


def find_implicants_qm(input_binstates, verbose=False):
    """ Finds the prime implicants (PI) using the Quine-McCluskey algorithm :cite:`Quine:1955`.

    Args:
        input_binstates (list / set) : A the binstates to condense.

    Returns:
        PI (set): a set of prime implicants.

    # Authors: Alex Gates
    """

    # we start with an empty set of implicants
    matched_implicants = set()
    done = False

    # repeat the following until no matches are found
    while not done:

        # split up the input_binstates into groups based on the number of 1s (density)
        density_groups = make_density_groups(input_binstates)

        # now clear everything for the new pass
        input_binstates = set()
        used = set()

        # Find the prime implicants

        # for each possible density
        for density in density_groups:
            # first make sure there are other binstates with the next density
            if density + 1 in density_groups:

                # then we pass through the binstates
                for binstate0 in density_groups[density]:

                    # An optimization due to Thomas Pircher, https://github.com/tpircher/quine-mccluskey/blob/master/quine_mccluskey/qm.py
                    # The Quine-McCluskey algorithm compares t1 with
                    # each element of the next group. (Normal approach)
                    # But in reality it is faster to construct all
                    # possible permutations of t1 by adding a '1' in
                    # opportune positions and check if this new term is
                    # contained in the set groups[key_next].

                    for idx, b0 in enumerate(binstate0):
                        if b0 == '0':
                            binstate1 = flip_binstate_bit(binstate0, idx)
                            if binstate1 in density_groups[density + 1]:
                                # keep track of the covered binary states
                                used.add(binstate0)
                                used.add(binstate1)
                                # keep the new wildcard binstate for the next round
                                input_binstates.add(replace_wildcard(binstate0, idx))

        # now add back the implicants that were not matched
        for groups in list(density_groups.values()):
            matched_implicants |= groups - used

        # finally, check if this pass actually compressed any terms
        # we finish when we cant make any further compressions
        if len(used) == 0:
            done = True

    # finish up by adding back all of the uncovered binary states
    prime_implicants = matched_implicants
    for groups in list(density_groups.values()):
        prime_implicants |= groups

    return prime_implicants


@cython.cfunc
cdef inline bint _is_wildcard_symbol(object symbol):
    return (
        symbol == WILDCARD_SYMBOL
        or symbol == SYMMETRIC_WILDCARD_SYMBOL
        or symbol == '2'
        or symbol == 2
    )


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef bint pi_covers_fast(object implicant, object binstate):
    """Fast coverage check for wildcard schemata."""
    cdef list imp_list
    cdef Py_ssize_t i, n
    cdef object symbol
    cdef str binstate_str = <str>binstate

    if isinstance(implicant, str):
        imp_list = list(implicant)
    else:
        imp_list = [x for x in implicant]

    n = len(imp_list)
    if len(binstate_str) != n:
        raise ValueError("Implicant and binstate must have the same length")

    for i in range(n):
        symbol = imp_list[i]
        if symbol == binstate_str[i]:
            continue
        if _is_wildcard_symbol(symbol):
            continue
        return False
    return True


@cython.cfunc
cdef list _normalize_index_groups(object groups):
    if not groups:
        return []
    cdef list normalized = []
    cdef object group
    cdef object idx
    for group in groups:
        if not group:
            normalized.append([])
            continue
        normalized.append([int(idx) for idx in group])
    return normalized


@cython.cfunc
cdef list _init_implicant_stack(object two_symbol):
    cdef list stack = []
    cdef object item
    if isinstance(two_symbol, str):
        stack.append(list(two_symbol))
    else:
        for item in two_symbol:
            if isinstance(item, str):
                stack.append(list(item))
            else:
                stack.append([c for c in item])
    return stack


@cython.cfunc
cdef void _permute_assign(
    list base_implicant,
    list idxs,
    list chars,
    Py_ssize_t pos,
    Py_ssize_t n,
    list stack,
    set seen,
):
    cdef Py_ssize_t i
    cdef object tmp
    cdef dict used

    if pos == n:
        new_implicant = base_implicant[:]
        for i in range(n):
            new_implicant[idxs[i]] = chars[i]
        tmp = tuple(new_implicant)
        if tmp not in seen:
            stack.append(new_implicant)
        return

    used = {}
    for i in range(pos, n):
        tmp = chars[i]
        if tmp in used:
            continue
        used[tmp] = None
        chars[pos], chars[i] = chars[i], chars[pos]
        _permute_assign(base_implicant, idxs, chars, pos + 1, n, stack, seen)
        chars[pos], chars[i] = chars[i], chars[pos]


@cython.cfunc
cdef void _enqueue_permutations(list base_implicant, list idxs, list stack, set seen):
    cdef Py_ssize_t n = len(idxs)
    if n <= 1:
        return
    chars = [base_implicant[i] for i in idxs]
    _permute_assign(base_implicant, idxs, chars, 0, n, stack, seen)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef list expand_ts_logic_fast(object two_symbols, object permut_indexes):
    """Generate all permutations for two-symbol schemata groups."""
    cdef list idx_groups = _normalize_index_groups(permut_indexes)
    if not idx_groups:
        return _init_implicant_stack(two_symbols)

    cdef list stack = _init_implicant_stack(two_symbols)
    cdef list results = []
    cdef set seen = set()
    cdef list implicant
    cdef tuple key
    cdef list idxs

    while stack:
        implicant = <list>stack.pop()
        key = tuple(implicant)
        if key in seen:
            continue
        seen.add(key)
        results.append(list(implicant))
        for idxs in idx_groups:
            if not idxs:
                continue
            _enqueue_permutations(implicant, idxs, stack, seen)

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef bint ts_covers_fast(object two_symbol, object permut_indexes, object binstate):
    """Fast coverage test for two-symbol schemata with permutations."""
    cdef list idx_groups = _normalize_index_groups(permut_indexes)
    if not idx_groups:
        return pi_covers_fast(two_symbol, binstate)

    cdef list stack = _init_implicant_stack(two_symbol)
    cdef set seen = set()
    cdef list implicant
    cdef tuple key
    cdef list idxs

    while stack:
        implicant = <list>stack.pop()
        key = tuple(implicant)
        if key in seen:
            continue
        seen.add(key)
        if pi_covers_fast(implicant, binstate):
            return True
        for idxs in idx_groups:
            if not idxs:
                continue
            _enqueue_permutations(implicant, idxs, stack, seen)

    return False


def __pi_covers(implicant, binstate):
    """Determine if a binarystate is covered by a specific implicant."""
    return pi_covers_fast(implicant, binstate)


def expand_wildcard_schemata(schemata):
    """
    Expand a wildcard schemata to list all binary states it covers.

    Args:
        schemata (string): the wildcard shemata
    Returns:
        binary_states (list): list of all binary states covered by the schemata
    """

    # count the number of wildcard symbols
    nwildcards = schemata.count(WILDCARD_SYMBOL)

    # if there arent any symbols, return the original schemata
    if nwildcards == 0:
        return [schemata]
    else:
        binary_states = []
        for wildstatenum in range(2**nwildcards):
            wildstates = statenum_to_binstate(wildstatenum, nwildcards)
            wnum = 0
            newstate = ''
            for b in schemata:
                if b == WILDCARD_SYMBOL:
                    newstate += wildstates[wnum]
                    wnum += 1
                else:
                    newstate += b
            binary_states.append(newstate)
        return binary_states


def return_pi_coverage(prime_implicants):
    """Computes the binary states coverage by Prime Implicant schematas.

    Args:
        prime_implicants (set): a set of prime implicants.
            This is returned by `find_implicants_qm`.
    Returns:
        pi_coverage (dict) : a dictionary of coverage where keys are input states and values are lists of the Prime Implicants covering that input.
    """

    pi_coverage = dict()
    for pi in prime_implicants:
        for binstate in expand_wildcard_schemata(pi):
            if binstate not in pi_coverage:
                pi_coverage[binstate] = set()
            pi_coverage[binstate].add(pi)

    return pi_coverage


def input_wildcard_coverage(pi_coverage):
    """Computes the binary states coverage by Prime Implicant schematas.

    Args:
        pi_coverage (dict): a dict mapping binary states to their prime implicants.
            This is returned by `return_pi_coverage`.
    Returns:
        input_wildcard_coverage (dict) : a dictionary of coverage where keys are inputs and values are lists of wither a WildCard covers that input.

    """
    # number of inputs
    k = len(next(iter(pi_coverage)))

    input_to_wildcards = {i: dict() for i in range(k)}
    for binstate, piset in pi_coverage.items():
        for i in range(k):
            input_to_wildcards[i][binstate] = tuple(pi[i] == WILDCARD_SYMBOL for pi in piset)

    return input_to_wildcards
