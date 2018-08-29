"""
Implementation of Algorithms 2 and 3 of [Ball2018]_

.. moduleauthor:: Fabian Ball <fabian.ball@kit.edu>
"""


def compute_orbit_partition(S, n):
    """
    Algorithm to compute the partition of orbits from a set of generators of the permutation group.

    Data representations are:
    * Permutation ``pi``: array-like, i.e. ``pi[i]`` corresponds to :math:`i \mapsto i^\pi`, which MUST be hashable
    * Set of generators ``S``: any container that is iterable repeatedly (e.g. ``list`` or ``set``)
    * Partition ``O``: array-like, i.e. ``O[i]`` corresponds to the (arbitrary) cluster id of node ``i``

    Besides the arbitrariness of cluster ids, this algorithm guarantees that ``O[i] == i`` if ``i`` is the smallest
    node id/label that is part of this cluster.

    Worst case time complexity is :math:`O(|S|\cdot n)`.

    Usage:

    >>> from sparsepermutation import SparsePermutation
    >>> from orbits import compute_orbit_partition
    >>> p1 = SparsePermutation([1, 0, 2, 3, 5, 4, 6, 7, 8, 9])  # (0 1)(4 5)
    >>> p2 = SparsePermutation([0, 1, 2, 3, 4, 5, 6, 7, 9, 8])  # (8 9)
    >>> p3 = SparsePermutation([0, 8, 2, 3, 4, 5, 6, 7, 1, 9])  # (1 8)
    >>> orbits = compute_orbit_partition([p1, p2, p3], 10)
    >>> assert orbits == [0, 0, 2, 3, 4, 4, 6, 7, 0, 0]

    :param S: A bunch of generators for a permutation group
    :type S: list | set | tuple
    :param n: The length of the set of elements the generated group acts on
    :type n: int
    :return: The partition of orbits of the permutation group
    :rtype: list
    """
    O = [-1] * n
    colored = 0
    U = set()
    for i in range(n):
        if O[i] != -1:
            continue

        N = set()
        colored += _color(O, i, S, i, N, U)
        if colored == n:
            return O

        while N:
            j = N.pop()
            colored += _color(O, j, S, i, N, U)
            if colored == n:
                return O

    return O


def _color(O, i, S, col, N, U):
    """
    A procedure that colors node ``i`` and all nodes that can be reached by permutations in ``S``.
    All newly reached nodes are added to ``N`` and every combination ``(i, pi)`` is added to ``U``.
    ``O`` is the partial orbit partition that is updated by coloring (= assigning orbit ids) the nodes.

    :param O: Partial orbit partition, i.e. not all nodes are assigned a cluster id, yet
    :type O: list
    :param i: Node id that shall be explored
    :type i: int
    :param S: Set of generators
    :type S: list | set | tuple
    :param col: The color that shall be used as node id
    :type col: int
    :param N: A set of node ids to which nodes on the same orbit as *i* are added
    :type N: set
    :param U: A set to keep track of already explored permutations for a given node id
    :return: The number of nodes that were colored
    :rtype: int
    """
    if O[i] == -1:
        O[i] = col
        colored = 1
    else:
        colored = 0

    for pi in S:
        if (i, pi) in U:
            continue

        U.add((i, pi))
        j = pi[i]

        while j != i:
            if O[j] == -1:
                O[j] = col
                N.add(j)
                colored += 1
                U.add((j, pi))
            else:
                assert O[j] == col

            j = pi[j]

    return colored
