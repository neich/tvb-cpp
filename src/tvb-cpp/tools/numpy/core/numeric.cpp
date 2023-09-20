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

#include "numeric.h"

using namespace tvb;

TArray1d convolve(const TArray1d& A, const TArray1d& B) {
    TArray1d a, b;
    if (A.size() < B.size()) {
        a = B;
        b = A;
    }
    else {
        a = A;
        b = B;
    }
    b.reverseInPlace();

    int len = b.size() + a.size() - 1;
    TArray1d result = TArray1d::Zero(len);
    int bs = b.size() - 1;
    for (int i = 0; i < len; ++i) {
        double v = 0.0;
        for (int j = 0; j < b.size(); ++j) {
            int a_index = j - bs + i;
            v += b(j) * (a_index >= 0 && a_index < a.size() ? a(a_index) : 0.0);
        }
        result[i] = v;
    }
    return result;

}

TArray1dc convolve(const TArray1dc& A, const TArray1dc& B) {
    TArray1dc a, b;
    if (A.size() < B.size()) {
        a = B;
        b = A;
    }
    else {
        a = A;
        b = B;
    }
    b.reverseInPlace();

    int len = b.size() + a.size() - 1;
    TArray1dc result = TArray1dc::Zero(len);
    int bs = b.size() - 1;
    for (int i = 0; i < len; ++i) {
        std::complex<double> v = 0.0;
        for (int j = 0; j < b.size(); ++j) {
            int a_index = j - bs + i;
            v += b(j) * (a_index >= 0 && a_index < a.size() ? a(a_index) : 0.0);
        }
        result[i] = v;
    }
    return result;
}


TArray1d poly(const TArray1d& seq_of_zeros) {
/*
    """
    Find the coefficients of a polynomial with the given sequence of roots.

    Returns the coefficients of the polynomial whose leading coefficient
    is one for the given sequence of zeros (multiple roots must be included
    in the sequence as many times as their multiplicity; see Examples).
    A square matrix (or array, which will be treated as a matrix) can also
    be given, in which case the coefficients of the characteristic polynomial
    of the matrix are returned.

    Parameters
    ----------
    seq_of_zeros : array_like, shape (N,) or (N, N)
        A sequence of polynomial roots, or a square array or matrix object.

    Returns
    -------
    c : ndarray
        1D array of polynomial coefficients from highest to lowest degree:

        ``c[0] * x**(N) + c[1] * x**(N-1) + ... + c[N-1] * x + c[N]``
        where c[0] always equals 1.

    Raises
    ------
    ValueError
        If input is the wrong shape (the input must be a 1-D or square
        2-D array).

    See Also
    --------
    polyval : Compute polynomial values.
    roots : Return the roots of a polynomial.
    polyfit : Least squares polynomial fit.
    poly1d : A one-dimensional polynomial class.

    Notes
    -----
    Specifying the roots of a polynomial still leaves one degree of
    freedom, typically represented by an undetermined leading
    coefficient. [1]_ In the case of this function, that coefficient -
    the first one in the returned array - is always taken as one. (If
    for some reason you have one other point, the only automatic way
    presently to leverage that information is to use ``polyfit``.)

    The characteristic polynomial, :math:`p_a(t)`, of an `n`-by-`n`
    matrix **A** is given by

        :math:`p_a(t) = \\mathrm{det}(t\\, \\mathbf{I} - \\mathbf{A})`,

    where **I** is the `n`-by-`n` identity matrix. [2]_

    References
    ----------
    .. [1] M. Sullivan and M. Sullivan, III, "Algebra and Trignometry,
       Enhanced With Graphing Utilities," Prentice-Hall, pg. 318, 1996.

    .. [2] G. Strang, "Linear Algebra and Its Applications, 2nd Edition,"
       Academic Press, pg. 182, 1980.

    Examples
    --------
    Given a sequence of a polynomial's zeros:

    >>> np.poly((0, 0, 0)) # Multiple root example
    array([1., 0., 0., 0.])

    The line above represents z**3 + 0*z**2 + 0*z + 0.

    >>> np.poly((-1./2, 0, 1./2))
    array([ 1.  ,  0.  , -0.25,  0.  ])

    The line above represents z**3 - z/4

    >>> np.poly((np.random.random(1)[0], 0, np.random.random(1)[0]))
    array([ 1.        , -0.77086955,  0.08618131,  0.        ]) # random

    Given a square array object:

    >>> P = np.array([[0, 1./3], [-1./2, 0]])
    >>> np.poly(P)
    array([1.        , 0.        , 0.16666667])

    Note how in all cases the leading coefficient is always 1.

    """
    seq_of_zeros = atleast_1d(seq_of_zeros)
    sh = seq_of_zeros.shape

    if len(sh) == 2 and sh[0] == sh[1] and sh[0] != 0:
        seq_of_zeros = eigvals(seq_of_zeros)
    elif len(sh) == 1:
        dt = seq_of_zeros.dtype
        # Let object arrays slip through, e.g. for arbitrary precision
        if dt != object:
            seq_of_zeros = seq_of_zeros.astype(mintypecode(dt.char))
    else:
        raise ValueError("input must be 1d or non-empty square 2d array.")
*/

    TArray1d a = TArray1d::Ones(1);
    if (seq_of_zeros.size() == 0)
        return a;
    // dt = seq_of_zeros.dtype
    for (unsigned k = 0; k < seq_of_zeros.size(); ++k) {
        TArray1d b(2);
        b[0] = 1.0;
        b[1] = -seq_of_zeros[k];
        a = convolve(a, b);
    }
//    if issubclass(a.dtype.type, NX.complexfloating):
//        //# if complex roots are all complex conjugates, the roots are real.
//        roots = NX.asarray(seq_of_zeros, complex)
//        if NX.all(NX.sort(roots) == NX.sort(roots.conjugate())):
//            a = a.real.copy()

    return a;
}

TArray1d poly(const TArray1dc& seq_of_zeros) {
/*
    """
    Find the coefficients of a polynomial with the given sequence of roots.

    Returns the coefficients of the polynomial whose leading coefficient
    is one for the given sequence of zeros (multiple roots must be included
    in the sequence as many times as their multiplicity; see Examples).
    A square matrix (or array, which will be treated as a matrix) can also
    be given, in which case the coefficients of the characteristic polynomial
    of the matrix are returned.

    Parameters
    ----------
    seq_of_zeros : array_like, shape (N,) or (N, N)
        A sequence of polynomial roots, or a square array or matrix object.

    Returns
    -------
    c : ndarray
        1D array of polynomial coefficients from highest to lowest degree:

        ``c[0] * x**(N) + c[1] * x**(N-1) + ... + c[N-1] * x + c[N]``
        where c[0] always equals 1.

    Raises
    ------
    ValueError
        If input is the wrong shape (the input must be a 1-D or square
        2-D array).

    See Also
    --------
    polyval : Compute polynomial values.
    roots : Return the roots of a polynomial.
    polyfit : Least squares polynomial fit.
    poly1d : A one-dimensional polynomial class.

    Notes
    -----
    Specifying the roots of a polynomial still leaves one degree of
    freedom, typically represented by an undetermined leading
    coefficient. [1]_ In the case of this function, that coefficient -
    the first one in the returned array - is always taken as one. (If
    for some reason you have one other point, the only automatic way
    presently to leverage that information is to use ``polyfit``.)

    The characteristic polynomial, :math:`p_a(t)`, of an `n`-by-`n`
    matrix **A** is given by

        :math:`p_a(t) = \\mathrm{det}(t\\, \\mathbf{I} - \\mathbf{A})`,

    where **I** is the `n`-by-`n` identity matrix. [2]_

    References
    ----------
    .. [1] M. Sullivan and M. Sullivan, III, "Algebra and Trignometry,
       Enhanced With Graphing Utilities," Prentice-Hall, pg. 318, 1996.

    .. [2] G. Strang, "Linear Algebra and Its Applications, 2nd Edition,"
       Academic Press, pg. 182, 1980.

    Examples
    --------
    Given a sequence of a polynomial's zeros:

    >>> np.poly((0, 0, 0)) # Multiple root example
    array([1., 0., 0., 0.])

    The line above represents z**3 + 0*z**2 + 0*z + 0.

    >>> np.poly((-1./2, 0, 1./2))
    array([ 1.  ,  0.  , -0.25,  0.  ])

    The line above represents z**3 - z/4

    >>> np.poly((np.random.random(1)[0], 0, np.random.random(1)[0]))
    array([ 1.        , -0.77086955,  0.08618131,  0.        ]) # random

    Given a square array object:

    >>> P = np.array([[0, 1./3], [-1./2, 0]])
    >>> np.poly(P)
    array([1.        , 0.        , 0.16666667])

    Note how in all cases the leading coefficient is always 1.

    """
    seq_of_zeros = atleast_1d(seq_of_zeros)
    sh = seq_of_zeros.shape

    if len(sh) == 2 and sh[0] == sh[1] and sh[0] != 0:
        seq_of_zeros = eigvals(seq_of_zeros)
    elif len(sh) == 1:
        dt = seq_of_zeros.dtype
        # Let object arrays slip through, e.g. for arbitrary precision
        if dt != object:
            seq_of_zeros = seq_of_zeros.astype(mintypecode(dt.char))
    else:
        raise ValueError("input must be 1d or non-empty square 2d array.")
*/

    TArray1dc a = TArray1dc::Ones(1);
    if (seq_of_zeros.size() > 0) {
        // dt = seq_of_zeros.dtype
        for (unsigned k = 0; k < seq_of_zeros.size(); ++k) {
            TArray1dc v(2);
            v[0] = 1.0;
            v[1] = -seq_of_zeros[k];
            a = convolve(a, v);
        }
    }

    TArray1d result(a.size());
    for (unsigned i = 0; i < a.size(); ++i)
        result(i) = a(i).real();

    return result;
}
