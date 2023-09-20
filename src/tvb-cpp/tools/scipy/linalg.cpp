//    Copyright 2020-2021 Ignacio Martín <ignacio.martin@udg.edu>
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//            http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#include "linalg.h"

using namespace tvb;

TMatrix companion(const TArray1d& a) {
/*
    """
    Create a companion matrix.

    Create the companion matrix [1]_ associated with the polynomial whose
    coefficients are given in `a`.

    Parameters
    ----------
    a : (N,) array_like
        1-D array of polynomial coefficients. The length of `a` must be
        at least two, and ``a[0]`` must not be zero.

    Returns
    -------
    c : (N-1, N-1) ndarray
        The first row of `c` is ``-a[1:]/a[0]``, and the first
        sub-diagonal is all ones.  The data-type of the array is the same
        as the data-type of ``1.0*a[0]``.

    Raises
    ------
    ValueError
        If any of the following are true: a) ``a.ndim != 1``;
        b) ``a.size < 2``; c) ``a[0] == 0``.

    Notes
    -----
    .. versionadded:: 0.8.0

    References
    ----------
    .. [1] R. A. Horn & C. R. Johnson, *TMatrix Analysis*.  Cambridge, UK:
        Cambridge University Press, 1999, pp. 146-7.

    Examples
    --------
    >>> from scipy.linalg import companion
    >>> companion([1, -10, 31, -30])
    array([[ 10., -31.,  30.],
           [  1.,   0.,   0.],
           [  0.,   1.,   0.]])

    """
*/
    // a = np.atleast_1d(a)

//    if a.ndim != 1:
//        raise ValueError("Incorrect shape for `a`.  `a` must be "
//                         "one-dimensional.")

    int n = a.size();
    if (a.size() < 2)
        throw("The length of `a` must be at least 2.");

    if (a[0] == 0)
        throw("The first coefficient in `a` must not be zero.");

    TVector first_row = -a(Eigen::seqN(1, n - 1)) / (1.0 * a[0]);
    TMatrix c = TMatrix::Zero(n - 1, n - 1); //np.zeros((n - 1, n - 1), dtype=first_row.dtype)
    c.row(0) = first_row;
    for (unsigned i = 0; i < n-2; ++i)
        c(i+1, i) = 1.0;
    // c[list(range(1, n - 1)), list(range(0, n - 2))] = 1
    return c;
}

