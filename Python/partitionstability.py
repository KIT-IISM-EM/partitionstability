"""
Implementation of Algorithms 1 and 4 of [Ball2018]_

.. moduleauthor:: Fabian Ball <fabian.ball@kit.edu>
"""
from orbits import compute_orbit_partition


def test_stability(P, S):
    """
    Test partition stability.

    Usage:

    >>> from sparsepermutation import SparsePermutation
    >>> from partitionstability import test_stability
    >>> p1 = SparsePermutation([1, 0, 2, 3, 5, 4, 6, 7, 8, 9])  # (0 1)(4 5)
    >>> p2 = SparsePermutation([0, 1, 2, 3, 4, 5, 6, 7, 9, 8])  # (8 9)
    >>> p3 = SparsePermutation([0, 8, 2, 3, 4, 5, 6, 7, 1, 9])  # (1 8)
    >>> S = [p1, p2, p3]
    >>> P = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    >>> assert test_stability(P, S) == True

    :param P: A partition in an array-like representation, clusters are identified by cluster ids
    :type P: list | tuple
    :param S: A set of generators for a permutation group
    :type S: list | tuple | set
    :return: True, if the partition is stable
    :rtype: bool
    """
    if not S:
        return True

    O = compute_orbit_partition(S, len(P))

    if geq(P, O):
        return True

    for pi in S:
        P_pi = [pi[i] for i in P]  # Apply pi on P
        if not geq(P, P_pi):
            return False

    return True


def geq(P, Q):
    """
    Test if ``P`` is coarser than or equal to ``Q``.
    Both, equality and coarseness are up to label isomorphism, i.e. partitions [0, 0, 1] and [3, 3, 7] are equal.
    A partition ``P`` is coarser than a partition ``Q`` if each cluster in ``Q`` is a subset of a cluster in ``P``.

    Data representation is:
     * Partition ``P``: array-like, i.e. ``P[i]`` corresponds to the (arbitrary) cluster id of node ``i``

    Usage:

    >>> from partitionstability import geq
    >>> P = [0,0,0,1,1,1]
    >>> P_prime = [1,1,1,0,0,0]
    >>> Q = [3,3,0,2,2,1]
    >>> R = [3,3,3,3,1,1]
    >>> assert geq(P, Q) == True
    >>> assert geq(P, P_prime) == True
    >>> assert geq(P, P_prime) == geq(P_prime, P)  # Cluster ids are arbitrary!
    >>> assert geq(P, R) == False
    >>> assert geq(R, P) == False
    >>> S = list(range(100))
    >>> from random import shuffle
    >>> shuffle(S)
    >>> assert geq(list(range(100)), S) == True

    :param P: A partition in an array-like representation of node ids
    :type P: list | tuple
    :param Q: A partition in an array-like representation of node ids
    :type Q: list | tuple
    :return: True, if P >= Q
    :rtype: bool
    """
    assert len(P) == len(Q)
    maps = {}  # Save cluster ids from Q in P. Multiple clusters in Q can map to the same cluster id in P

    for i in range(len(P)):
        if Q[i] in maps:  # The orbit id in Q was already seen
            o_id = maps[Q[i]]
        else:  # The first occurrence of an orbit id is saved
            o_id = P[i]  # The orbit id in P
            maps[Q[i]] = o_id  # The orbit id in Q maps to the one in P

        if P[i] != o_id:  # The orbit ids must match
            return False

    return True
