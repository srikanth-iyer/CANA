# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False
"""Cython-accelerated helpers for ``cana.boolean_node``."""

import cython

from cpython.dict cimport PyDict_Size
from cpython.list cimport (
    PyList_Check,
    PyList_GET_ITEM,
    PyList_GET_SIZE,
)
from cpython.tuple cimport PyTuple_Check, PyTuple_GET_ITEM
from cpython.sequence cimport PySequence_Size
from cpython.exc cimport PyErr_Clear


@cython.cfunc
cdef inline Py_ssize_t _safe_sequence_size(object obj) except -1:
    cdef Py_ssize_t size = PySequence_Size(obj)
    if size < 0:
        PyErr_Clear()
        return 0
    return size


@cython.cfunc
cdef inline double _sum_group_lengths(list f_theta) except? -1.0:
    cdef Py_ssize_t len_f_theta = PyList_GET_SIZE(f_theta)
    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef Py_ssize_t group_len
    cdef double total = 0.0
    cdef object ts_obj
    cdef object groups_obj
    cdef object group_obj
    cdef tuple ts
    cdef list groups

    if len_f_theta <= 0:
        return 0.0

    for i in range(len_f_theta):
        ts_obj = <object>PyList_GET_ITEM(f_theta, i)
        if PyTuple_Check(ts_obj):
            ts = <tuple>ts_obj
        else:
            ts = tuple(ts_obj)
        groups_obj = <object>PyTuple_GET_ITEM(ts, 1)
        if groups_obj is None or not groups_obj:
            continue
        if not PyList_Check(groups_obj):
            groups_obj = list(groups_obj)
        groups = <list>groups_obj
        for j in range(PyList_GET_SIZE(groups)):
            group_obj = <object>PyList_GET_ITEM(groups, j)
            if PyList_Check(group_obj):
                group_len = PyList_GET_SIZE(group_obj)
            else:
                group_len = _safe_sequence_size(group_obj)
            total += group_len

    return total


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef double input_symmetry_mean_fast(dict ts_coverage, int k):
    """Fast computation of :meth:`BooleanNode.input_symmetry_mean`.

    Parameters
    ----------
    ts_coverage: dict
        Mapping of LUT entries to the list of covering Two-Symbol schemata.
    k: int
        Input degree of the Boolean node.

    Returns
    -------
    float
        Mean number of permutable inputs over the LUT entries.
    """
    cdef double numerator = 0.0
    cdef double inv_total_states
    cdef double inv_len
    cdef Py_ssize_t len_f_theta
    cdef list f_theta

    if ts_coverage is None:
        return 0.0
    if PyDict_Size(ts_coverage) == 0:
        return 0.0
    if k <= 0:
        return 0.0

    inv_total_states = 1.0 / (1 << k)

    for f_theta in ts_coverage.values():
        if not f_theta:
            continue
        if not PyList_Check(f_theta):
            f_theta = list(f_theta)
        f_theta = <list>f_theta
        len_f_theta = PyList_GET_SIZE(f_theta)
        if len_f_theta == 0:
            continue
        inv_len = 1.0 / len_f_theta
        numerator += _sum_group_lengths(f_theta) * inv_len

    return numerator * inv_total_states


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef double input_symmetry_mean_annigen_fast(dict ts_coverage):
    """Fast computation of :meth:`BooleanNode.input_symmetry_mean_anni_gen`."""
    cdef double numerator = 0.0
    cdef double inv_len
    cdef Py_ssize_t len_f_theta
    cdef Py_ssize_t coverage_count = 0
    cdef list f_theta

    if ts_coverage is None:
        return 0.0
    if PyDict_Size(ts_coverage) == 0:
        return 0.0

    for f_theta in ts_coverage.values():
        if not f_theta:
            continue
        if not PyList_Check(f_theta):
            f_theta = list(f_theta)
        f_theta = <list>f_theta
        len_f_theta = PyList_GET_SIZE(f_theta)
        if len_f_theta == 0:
            continue
        inv_len = 1.0 / len_f_theta
        numerator += _sum_group_lengths(f_theta) * inv_len
        coverage_count += 1

    if coverage_count == 0:
        return 0.0
    return numerator / coverage_count
