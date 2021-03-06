"""
Simple implementation of a sparse permutation datastructure which is more efficient than
saving permutation lists of length ``n`` each.

.. moduleauthor:: Fabian Ball <fabian.ball@kit.edu>
"""

__all__ = ['SparsePermutation']


class SparsePermutationArray(dict):
    """
    Simple implementation of a dict which returns 'identity' if a key is missing.
    This is used to represent a permutation in a sparse way (cycle notation).
    However, no checks are performed in terms of validity of the cycles.

    For internal use only, use ``SparsePermutation`` to create a sparse permutation
    """
    def __missing__(self, key):
        return key

    def __setitem__(self, key, value):
        if key == value:
            if key in self:
                self.__delattr__(key)
        else:
            super(SparsePermutationArray, self).__setitem__(key, value)

    def __str__(self):
        # Produce a nice cycle representation
        keys = set(self.keys())
        ret = ''

        while keys:
            k = min(keys)
            keys.remove(k)
            cycle = [k]
            img_k = self[k]
            while img_k != k:
                cycle.append(img_k)
                keys.remove(img_k)
                img_k = self[img_k]

            ret += '({})'.format(' '.join(map(str, cycle)))

        return ret


class SparsePermutation(object):
    """
    A minimal implementation of a sparse permutation that partly mimics the Permutation
    implementation of ``sympy.combinatorics.Permutation``.

    Usage:

    >>> from sparsepermutation import SparsePermutation
    >>> p = [0, 1, 3, 2]  # Permutation (2 3)
    >>> sp = SparsePermutation(p)  # Sparsify
    >>> q = list(sp)  # Convert back to list
    >>> assert p == q
    """
    def __init__(self, permutation_array):
        """
        Create a sparse permutation given a permutation list of length ``n``.

        :param permutation_array: Array-like permutation representation
        :type permutation_array: list | tuple
        """
        self._perm = SparsePermutationArray()
        self._n = len(permutation_array)  # Save the length of the permutation
        self._support = 0
        for k, v in enumerate(permutation_array):
            if k != v:
                self._perm[k] = v
                self._support += 1  # Increase the support count of the permutation

    def __call__(self, *args):
        """
        Convenient application of this permutation.
        """
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, int):  # Return the image
                return self[arg]
            else:  # Assume an iterable that is indexable, e.g. a list
                try:
                    t = type(arg)
                    return t([arg[self[i]] for i in range(len(arg))])
                except TypeError:
                    raise TypeError('The argument must be iterable and indexable')
        elif len(args) > 1:  # Apply permutation on every argument
            return tuple(self(arg) for arg in args)
        else:
            raise TypeError('No argument given')

    def __getitem__(self, item):
        if item >= self._n:
            raise KeyError(item)
        else:
            return self._perm[item]

    def __iter__(self):
        for i in range(self._n):
            yield self._perm[i]

    def __hash__(self):
        return hash(tuple(self))

    def __str__(self):
        return str(self._perm)

    #############
    # sympy.combinatorics.Permutation mimics
    #############
    def length(self):
        """
        :return: The number of affected elements
        :rtype: int
        """
        return self._support

    def atoms(self):
        """
        :return: The set of elements
        :rtype: set
        """
        return set(range(self._n))

    def support(self):
        """
        :return: The list of affected elements
        :rtype: list
        """
        return list(self._perm)

    @property
    def size(self):
        """
        :return: The number of elements the permutation acts on
        :rtype: int
        """
        return self._n
