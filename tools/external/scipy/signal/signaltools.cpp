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

#include <unordered_set>

#include <definitions.h>

#include <external/numpy/numpy.h>
#include <external/scipy/linalg.h>

#include <external/scipy/signal/signaltools.h>

using namespace tvb;
using namespace std;

Vectord detrend_linear(const tvb::Vectord& data) {
    double xmean, ymean;
    double temp;
    double Sxy;
    double Sxx;

    double grad;
    double yint;

    int m = data.size();
    Vectord x(m);
    Vectord y = data;

    /********************************
    Set the X axis Liner Values
    *********************************/
    for (unsigned i = 0; i < m; i++)
        x[i] = i;

    /********************************
    Calculate the mean of x and y
    *********************************/
    xmean = 0;
    ymean = 0;
    for (unsigned i = 0; i < m; i++)
    {
        xmean += x[i];
        ymean += y[i];
    }
    xmean /= m;
    ymean /= m;

    /********************************
    Calculate Covariance
    *********************************/
    temp = 0;
    for (unsigned i = 0; i < m; i++)
        temp += x[i] * y[i];
    Sxy = temp / m - xmean * ymean;

    temp = 0;
    for (unsigned i = 0; i < m; i++)
        temp += x[i] * x[i];
    Sxx = temp / m - xmean * xmean;

    /********************************
    Calculate Gradient and Y intercept
    *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    /********************************
    Removing Linear Trend
    *********************************/
    for (unsigned i = 0; i < m; i++)
        y[i] = y[i] - (grad * i + yint);

    return y;
}


//Vectord detrend(const Vectord& data, const string& type, const Vectori& bp=, bool overwrite_data) {
////    """
////    Remove linear trend along axis from data.
////
////    Parameters
////    ----------
////    data : array_like
////        The input data.
////    axis : int, optional
////        The axis along which to detrend the data. By default this is the
////        last axis (-1).
////    type : {"linear", "constant"}, optional
////        The type of detrending. If ``type == "linear"`` (default),
////        the result of a linear least-squares fit to `data` is subtracted
////        from `data`.
////        If ``type == "constant"``, only the mean of `data` is subtracted.
////    bp : array_like of ints, optional
////        A sequence of break points. If given, an individual linear fit is
////        performed for each part of `data` between two break points.
////        Break points are specified as indices into `data`. This parameter
////        only has an effect when ``type == "linear"``.
////    overwrite_data : bool, optional
////        If True, perform in place detrending and avoid a copy. Default is False
////
////    Returns
////    -------
////    ret : ndarray
////        The detrended input data.
////
////    Examples
////    --------
////    >>> from scipy import signal
////    >>> randgen = np.random.RandomState(9)
////    >>> npoints = 1000
////    >>> noise = randgen.randn(npoints)
////    >>> x = 3 + 2*np.linspace(0, 1, npoints) + noise
////    >>> (signal.detrend(x) - noise).max() < 0.01
////    True
////
////    """
//
//    static const unordered_set<string> types = {"linear", "l", "constant", "c"};
//
//    if (types.count(type) != 1)
//        throw("Trend type must be \"linear\" or \"constant\".");
//    // data = np.asarray(data)
////    dtype = data.dtype.char
////    if dtype not in "dfDF":
////        dtype = "d"
//    if (type == "constant" || type == "c") {
//        Vectord ret = data - data.sum()/data.size();
//        return ret;
//    }
//    else {
////        dshape = data.shape
//        int N = data.size();
//        auto bps = bp;
//        sort(bps.begin(), bps.end());
////        bp = np.sort(np.unique(np.r_[0, bp, N]))
//        if ((bps > N).any())
//            throw("Breakpoints must be less than length of data along given axis.");
//        int Nreg = bps.size() - 1;
////        # Restructure data so that axis is along first dimension and
////        #  all other dimensions are collapsed into second dimension
//        int rnk = 1;
//        if axis < 0:
//            axis = axis + rnk
//        newdims = np.r_[axis, 0:axis, axis + 1:rnk]
//        newdata = np.reshape(np.transpose(data, tuple(newdims)),
//                             (N, _prod(dshape) // N))
//        if not overwrite_data:
//            newdata = newdata.copy()  # make sure we have a copy
//        if newdata.dtype.char not in "dfDF":
//            newdata = newdata.astype(dtype)
//        # Find leastsq fit and remove it for each piece
//        for m in range(Nreg):
//            Npts = bp[m + 1] - bp[m]
//            A = np.ones((Npts, 2), dtype)
//            A[:, 0] = np.cast[dtype](np.arange(1, Npts + 1) * 1.0 / Npts)
//            sl = slice(bp[m], bp[m + 1])
//            coef, resids, rank, s = linalg.lstsq(A, newdata[sl])
//            newdata[sl] = newdata[sl] - np.dot(A, coef)
//        # Put data back in original shape.
//        tdshape = np.take(dshape, newdims, 0)
//        ret = np.reshape(newdata, tuple(tdshape))
//        vals = list(range(1, rnk))
//        olddims = vals[:axis] + [0] + vals[axis:]
//        ret = np.transpose(ret, tuple(olddims))
//        return ret


Vectord even_ext(const Vectord&x, int n) {
/*
    """
    Even extension at the boundaries of an array

    Generate a new ndarray by making an even extension of `x` along an axis.

    Parameters
    ----------
    x : ndarray
        The array to be extended.
    n : int
        The number of elements by which to extend `x` at each end of the axis.
    axis : int, optional
        The axis along which to extend `x`. Default is -1.

    Examples
    --------
    >>> from scipy.signal._arraytools import even_ext
    >>> a = np.array([[1, 2, 3, 4, 5], [0, 1, 4, 9, 16]])
    >>> even_ext(a, 2)
    array([[ 3,  2,  1,  2,  3,  4,  5,  4,  3],
           [ 4,  1,  0,  1,  4,  9, 16,  9,  4]])

    Even extension is a "mirror image" at the boundaries of the original array:

    >>> t = np.linspace(0, 1.5, 100)
    >>> a = 0.9 * np.sin(2 * np.pi * t**2)
    >>> b = even_ext(a, 40)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(arange(-40, 140), b, 'b', lw=1, label='even extension')
    >>> plt.plot(arange(100), a, 'r', lw=2, label='original')
    >>> plt.legend(loc='best')
    >>> plt.show()
    """
*/
    if (n < 1)
        return x;
    if (n > x.size() - 1)
        throw (string_format("The extension length n (%d) is too big. It must not exceed x.shape[axis]-1, which is %d.", n,
                      x.size() - 1));
    Vectord left_ext = x(Eigen::seqN(n, n, -1));
    // left_ext = axis_slice(x, start=n, stop=0, step=-1, axis=axis)
    Vectord right_ext = x(Eigen::seqN(Eigen::last, n, -1));
    // right_ext = axis_slice(x, start=-2, stop=-(n + 2), step=-1, axis=axis)
//    ext = np.concatenate((left_ext,
//                          x,
//                          right_ext),
//                         axis=axis)
    Vectord ext(x.size() + 2 * n);
    ext << left_ext, x, right_ext;
    return ext;
}


Vectord const_ext(const Vectord& x, int n) {
/*
    """
    Constant extension at the boundaries of an array

    Generate a new ndarray that is a constant extension of `x` along an axis.

    The extension repeats the values at the first and last element of
    the axis.

    Parameters
    ----------
    x : ndarray
        The array to be extended.
    n : int
        The number of elements by which to extend `x` at each end of the axis.
    axis : int, optional
        The axis along which to extend `x`. Default is -1.

    Examples
    --------
    >>> from scipy.signal._arraytools import const_ext
    >>> a = np.array([[1, 2, 3, 4, 5], [0, 1, 4, 9, 16]])
    >>> const_ext(a, 2)
    array([[ 1,  1,  1,  2,  3,  4,  5,  5,  5],
           [ 0,  0,  0,  1,  4,  9, 16, 16, 16]])

    Constant extension continues with the same values as the endpoints of the
    array:

    >>> t = np.linspace(0, 1.5, 100)
    >>> a = 0.9 * np.sin(2 * np.pi * t**2)
    >>> b = const_ext(a, 40)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(arange(-40, 140), b, 'b', lw=1, label='constant extension')
    >>> plt.plot(arange(100), a, 'r', lw=2, label='original')
    >>> plt.legend(loc='best')
    >>> plt.show()
    """
*/
    if (n < 1)
        return x;
    Vectord left_ext = Vectord::Constant(x(0), n);
    Vectord right_ext = Vectord::Constant(x(Eigen::last), n);
    Vectord ext(x.size() + 2 * n);
    ext << left_ext, x, right_ext;
    return ext;
}

Vectord odd_ext(const Vectord&x, int n) {
    /*  """
      Odd extension at the boundaries of an array

      Generate a new ndarray by making an odd extension of `x` along an axis.

      Parameters
      ----------
      x : ndarray
          The array to be extended.
      n : int
          The number of elements by which to extend `x` at each end of the axis.
      axis : int, optional
          The axis along which to extend `x`. Default is -1.

      Examples
      --------
      >>> from scipy.signal._arraytools import odd_ext
      >>> a = np.array([[1, 2, 3, 4, 5], [0, 1, 4, 9, 16]])
      >>> odd_ext(a, 2)
      array([[-1,  0,  1,  2,  3,  4,  5,  6,  7],
             [-4, -1,  0,  1,  4,  9, 16, 23, 28]])

      Odd extension is a "180 degree rotation" at the endpoints of the original
      array:

      >>> t = np.linspace(0, 1.5, 100)
      >>> a = 0.9 * np.sin(2 * np.pi * t**2)
      >>> b = odd_ext(a, 40)
      >>> import matplotlib.pyplot as plt
      >>> plt.plot(arange(-40, 140), b, 'b', lw=1, label='odd extension')
      >>> plt.plot(arange(100), a, 'r', lw=2, label='original')
      >>> plt.legend(loc='best')
      >>> plt.show()
      """*/
    if (n < 1)
        return x;
    if (n > x.size() - 1)
        throw (string_format("The extension length n (%d) is too big. It must not exceed x.shape[axis]-1, which is %d.",
                      n, x.size() - 1));
    double left_end = x[0]; // axis_slice(x, start=0, stop=1, axis=axis)
    Vectord left_ext = x(Eigen::seqN(n, n, -1)); // axis_slice(x, start=n, stop=0, step=-1, axis=axis)
    double right_end = x(Eigen::last); // axis_slice(x, start=-1, axis=axis)
    Vectord right_ext = x(Eigen::seqN(Eigen::last - 1, n, -1)); // axis_slice(x, start=-2, stop=-(n + 2), step=-1, axis=axis)
//    ext = np.concatenate((2 * left_end - left_ext,
//                                 x,
//                                 2 * right_end - right_ext),
//                         axis = axis)
    Vectord ext(x.size() + 2 * n);
    ext(Eigen::seqN(0, n)) = 2*left_end - left_ext;
    ext(Eigen::seqN(n, x.size())) = x;
    ext(Eigen::seqN(n+x.size(), n)) = 2*right_end - right_ext;
    // ext << 2*left_end - left_ext, x, 2*right_end - right_ext;
    return ext;
}

pair<int, Vectord> _validate_pad(const string& padtype, int padlen, const Vectord& x, int ntaps) {
    // """Helper to validate padding for filtfilt"""
    static const unordered_set<string> pad_types = {"even", "odd", "constant", "none"};
    if (pad_types.count(padtype) != 1)
        throw ("Padtype must be 'even', 'odd', 'constant', or none.");

    if (padtype == "none")
        padlen = 0;

    int edge = padlen;

    // # x's 'axis' dimension must be bigger than edge.
    if (x.size() <= edge)
        throw (string_format("The length of the input vector x must be greater than padlen, which is %d.", edge));

    Vectord ext;
    if (padtype != "none" && edge > 0) {
//        # Make an extension of length `edge` at each
//        # end of the input array.
        if (padtype == "even")
            ext = even_ext(x, edge);
        else if (padtype == "odd")
            ext = odd_ext(x, edge);
        else
            ext = const_ext(x, edge);
    } else
        ext = x;
    return {edge, ext};
}

Vectord lfilter_zi(const Vectord& B, const Vectord& A) {
    /*  """
      Construct initial conditions for lfilter for step response steady-state.

      Compute an initial state `zi` for the `lfilter` function that corresponds
      to the steady state of the step response.

      A typical use of this function is to set the initial state so that the
      output of the filter starts at the same value as the first element of
      the signal to be filtered.

      Parameters
      ----------
      b, a : array_like (1-D)
          The IIR filter coefficients. See `lfilter` for more
          information.

      Returns
      -------
      zi : 1-D ndarray
          The initial state for the filter.

      See Also
      --------
      lfilter, lfiltic, filtfilt

      Notes
      -----
      A linear filter with order m has a state space representation (A, B, C, D),
      for which the output y of the filter can be expressed as::

          z(n+1) = A*z(n) + B*x(n)
          y(n)   = C*z(n) + D*x(n)

      where z(n) is a vector of length m, A has shape (m, m), B has shape
      (m, 1), C has shape (1, m) and D has shape (1, 1) (assuming x(n) is
      a scalar).  lfilter_zi solves::

          zi = A*zi + B

      In other words, it finds the initial condition for which the response
      to an input of all ones is a constant.

      Given the filter coefficients `a` and `b`, the state space matrices
      for the transposed direct form II implementation of the linear filter,
      which is the implementation used by scipy.signal.lfilter, are::

          A = scipy.linalg.companion(a).T
          B = b[1:] - a[1:]*b[0]

      assuming `a[0]` is 1.0; if `a[0]` is not 1, `a` and `b` are first
      divided by a[0].

      Examples
      --------
      The following code creates a lowpass Butterworth filter. Then it
      applies that filter to an array whose values are all 1.0; the
      output is also all 1.0, as expected for a lowpass filter.  If the
      `zi` argument of `lfilter` had not been given, the output would have
      shown the transient signal.

      >>> from numpy import array, ones
      >>> from scipy.signal import lfilter, lfilter_zi, butter
      >>> b, a = butter(5, 0.25)
      >>> zi = lfilter_zi(b, a)
      >>> y, zo = lfilter(b, a, ones(10), zi=zi)
      >>> y
      array([1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

      Another example:

      >>> x = array([0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0])
      >>> y, zf = lfilter(b, a, x, zi=zi*x[0])
      >>> y
      array([ 0.5       ,  0.5       ,  0.5       ,  0.49836039,  0.48610528,
          0.44399389,  0.35505241])

      Note that the `zi` argument to `lfilter` was computed using
      `lfilter_zi` and scaled by `x[0]`.  Then the output `y` has no
      transient until the input drops from 0.5 to 0.0.

      """

      # FIXME: Can this function be replaced with an appropriate
      # use of lfiltic?  For example, when b,a = butter(N,Wn),
      #    lfiltic(b, a, y=numpy.ones_like(a), x=numpy.ones_like(b)).
      #

      # We could use scipy.signal.normalize, but it uses warnings in
      # cases where a ValueError is more appropriate, and it allows
      # b to be 2D.*/

    // b = np.atleast_1d(b)
//    if b.ndim != 1:
//        raise ValueError("Numerator b must be 1-D.")
//    a = np.atleast_1d(a)
//    if a.ndim != 1:
//        raise ValueError("Denominator a must be 1-D.")

    auto a = A;
    auto b = B;
    while (a.size() > 1 && a[0] == 0.0)
        a = a(Eigen::seq(1, Eigen::last)); // [1:]
    if (a.size() < 1)
        throw ("There must be at least one nonzero `a` coefficient.");

    if (a[0] != 1.0) {
        // # Normalize the coefficients so a[0] == 1.
        b = b / a[0];
        a = a / a[0];
    }

    auto n = max(a.size(), b.size());

    // # Pad a or b with zeros so they are the same length.
    if (a.size() < n)
        a.conservativeResize(n);
    else if (b.size() < n)
        b.conservativeResize(n);

    Eigen::MatrixXd IminusA = Eigen::MatrixXd::Identity(n - 1, n - 1) - companion(a).transpose();
    Eigen::VectorXd bb = b(Eigen::seqN(1, Eigen::last)) - a(Eigen::seqN(1, Eigen::last)) * b[0];
    // # Solve zi = A*zi + B
    Vectord zi = IminusA.colPivHouseholderQr().solve(bb).array();
    // zi = np.linalg.solve(IminusA, B)

//    # For future reference: we could also use the following
//    # explicit formulas to solve the linear system:
//    #
//    # zi = np.zeros(n - 1)
//    # zi[0] = B.sum() / IminusA[:,0].sum()
//    # asum = 1.0
//    # csum = 0.0
//    # for k in range(1,n-1):
//    #     asum += a[k]
//    #     csum += b[k] - a[k]*b[0]
//    #     zi[k] = asum*zi[0] - csum

    return zi;
}


Vectord _linear_filter(const Vectord &B, const Vectord &A, const Vectord &x, const Vectord& zi={})
{
    /* normalize the filter coefs only once. */
    Vectord b = B / A(0);
    Vectord a = A / A(0);
    Vectord y(x.size());
    Vectord z = zi;

    for (unsigned k = 0; k < x.size(); k++) {
        if (b.size() > 1) {
            y[k] = z[0] + b[0] * x[k];
            /* Fill in middle delays */
            for (unsigned n = 1; n < b.size() - 1; n++) {
                z[n-1] = z[n] + x[k] * b[n] - y[k] * a[n];
            }
            /* Calculate last delay */
            z[z.size()-1] = x[k] * b(Eigen::last) - y[k] * a(Eigen::last);
        } else {
            y[k] = x[k] * b[0];
        }
    }
    return y;
}


Vectord lfilter(const Vectord& b, const Vectord& a, const Vectord& x, const Vectord& zi={}) {
/*
    """
    Filter data along one-dimension with an IIR or FIR filter.

    Filter a data sequence, `x`, using a digital filter.  This works for many
    fundamental data types (including Object type).  The filter is a direct
    form II transposed implementation of the standard difference equation
    (see Notes).

    The function `sosfilt` (and filter design using ``output='sos'``) should be
    preferred over `lfilter` for most filtering tasks, as second-order sections
    have fewer numerical problems.

    Parameters
    ----------
    b : array_like
        The numerator coefficient vector in a 1-D sequence.
    a : array_like
        The denominator coefficient vector in a 1-D sequence.  If ``a[0]``
        is not 1, then both `a` and `b` are normalized by ``a[0]``.
    x : array_like
        An N-dimensional input array.
    axis : int, optional
        The axis of the input data array along which to apply the
        linear filter. The filter is applied to each subarray along
        this axis.  Default is -1.
    zi : array_like, optional
        Initial conditions for the filter delays.  It is a vector
        (or array of vectors for an N-dimensional input) of length
        ``max(len(a), len(b)) - 1``.  If `zi` is None or is not given then
        initial rest is assumed.  See `lfiltic` for more information.

    Returns
    -------
    y : array
        The output of the digital filter.
    zf : array, optional
        If `zi` is None, this is not returned, otherwise, `zf` holds the
        final filter delay values.

    See Also
    --------
    lfiltic : Construct initial conditions for `lfilter`.
    lfilter_zi : Compute initial state (steady state of step response) for
                 `lfilter`.
    filtfilt : A forward-backward filter, to obtain a filter with linear phase.
    savgol_filter : A Savitzky-Golay filter.
    sosfilt: Filter data using cascaded second-order sections.
    sosfiltfilt: A forward-backward filter using second-order sections.

    Notes
    -----
    The filter function is implemented as a direct II transposed structure.
    This means that the filter implements::

       a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
                             - a[1]*y[n-1] - ... - a[N]*y[n-N]

    where `M` is the degree of the numerator, `N` is the degree of the
    denominator, and `n` is the sample number.  It is implemented using
    the following difference equations (assuming M = N)::

         a[0]*y[n] = b[0] * x[n]               + d[0][n-1]
           d[0][n] = b[1] * x[n] - a[1] * y[n] + d[1][n-1]
           d[1][n] = b[2] * x[n] - a[2] * y[n] + d[2][n-1]
         ...
         d[N-2][n] = b[N-1]*x[n] - a[N-1]*y[n] + d[N-1][n-1]
         d[N-1][n] = b[N] * x[n] - a[N] * y[n]

    where `d` are the state variables.

    The rational transfer function describing this filter in the
    z-transform domain is::

                             -1              -M
                 b[0] + b[1]z  + ... + b[M] z
         Y(z) = -------------------------------- X(z)
                             -1              -N
                 a[0] + a[1]z  + ... + a[N] z

    Examples
    --------
    Generate a noisy signal to be filtered:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    >>> t = np.linspace(-1, 1, 201)
    >>> x = (np.sin(2*np.pi*0.75*t*(1-t) + 2.1) +
    ...      0.1*np.sin(2*np.pi*1.25*t + 1) +
    ...      0.18*np.cos(2*np.pi*3.85*t))
    >>> xn = x + np.random.randn(len(t)) * 0.08

    Create an order 3 lowpass butterworth filter:

    >>> b, a = signal.butter(3, 0.05)

    Apply the filter to xn.  Use lfilter_zi to choose the initial condition of
    the filter:

    >>> zi = signal.lfilter_zi(b, a)
    >>> z, _ = signal.lfilter(b, a, xn, zi=zi*xn[0])

    Apply the filter again, to have a result filtered at an order the same as
    filtfilt:

    >>> z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])

    Use filtfilt to apply the filter:

    >>> y = signal.filtfilt(b, a, xn)

    Plot the original signal and the various filtered versions:

    >>> plt.figure
    >>> plt.plot(t, xn, 'b', alpha=0.75)
    >>> plt.plot(t, z, 'r--', t, z2, 'r', t, y, 'k')
    >>> plt.legend(('noisy signal', 'lfilter, once', 'lfilter, twice',
    ...             'filtfilt'), loc='best')
    >>> plt.grid(True)
    >>> plt.show()

    """
    a = np.atleast_1d(a)
*/
    if (a.size() == 1) {
        // TODO: implement this branch
//        # This path only supports types fdgFDGO to mirror _linear_filter below.
//        # Any of b, a, x, or zi can set the dtype, but there is no default
//        # casting of other types; instead a NotImplementedError is raised.
//        b = np.asarray(b)
//        a = np.asarray(a)
//        if b.ndim != 1 and a.ndim != 1:
//            raise ValueError('object of too small depth for desired array')
//        x = _validate_x(x)
//        inputs = [b, a, x]
//        if (zi.size() > 0) {
//            # _linear_filter does not broadcast zi, but does do expansion of
//            # singleton dims.
//            zi = np.asarray(zi)
//            if zi.ndim != x.ndim:
//                raise ValueError('object of too small depth for desired array')
//            expected_shape = list(x.shape)
//            expected_shape[axis] = b.shape[0] - 1
//            expected_shape = tuple(expected_shape)
////            # check the trivial case where zi is the right shape first
//            if zi.shape != expected_shape:
//                strides = zi.ndim * [None]
//                if axis < 0:
//                    axis += zi.ndim
//                for k in range(zi.ndim):
//                    if k == axis and zi.shape[k] == expected_shape[k]:
//                        strides[k] = zi.strides[k]
//                    elif k != axis and zi.shape[k] == expected_shape[k]:
//                        strides[k] = zi.strides[k]
//                    elif k != axis and zi.shape[k] == 1:
//                        strides[k] = 0
//                    else:
//                        raise ValueError('Unexpected shape for zi: expected '
//                                         '%s, found %s.' %
//                                         (expected_shape, zi.shape))
//                zi = np.lib.stride_tricks.as_strided(zi, expected_shape,
//                                                     strides)
//            inputs.append(zi)
//        dtype = np.result_type(*inputs)
//
//        if dtype.char not in 'fdgFDGO':
//            raise NotImplementedError("input type '%s' not supported" % dtype)

//        b = np.array(b, dtype=dtype)
//        a = np.array(a, dtype=dtype, copy=False)
//        auto bb = b / a[0];
//        x = np.array(x, dtype=dtype, copy=False)

//        auto out_full = convolve(b, x);
//        ind = out_full.ndim * [slice(None)];
//        if (zi.size() > 0) {
//            ind[axis] = slice(zi.shape[axis]);
//            out_full[tuple(ind)] += zi;
//        }
//
//        ind[axis] = slice(out_full.shape[axis] - len(b) + 1);
//        out = out_full[tuple(ind)];
//
//        if (zi.size() == 0)
//            return out;
//        else {
//            ind[axis] = slice(out_full.shape[axis] - len(b) + 1, None);
//            zf = out_full[tuple(ind)];
//            return out, zf;
//        }
    } else {
        if (zi.size() == 0)
            return _linear_filter(b, a, x);
        else
            // TODO: Implement initial state
            return _linear_filter(b, a, x, zi);
    }
    throw std::runtime_error("Not implemented");
}


Vectord filtfilt_pad(const Vectord& b, const Vectord& a, const Vectord& x,
         int padlen, const string& padtype) {
   /* """
    Apply a digital filter forward and backward to a signal.

    This function applies a linear digital filter twice, once forward and
    once backwards.  The combined filter has zero phase and a filter order
    twice that of the original.

    The function provides options for handling the edges of the signal.

    The function `sosfiltfilt` (and filter design using ``output='sos'``)
    should be preferred over `filtfilt` for most filtering tasks, as
    second-order sections have fewer numerical problems.

    Parameters
    ----------
    b : (N,) array_like
        The numerator coefficient vector of the filter.
    a : (N,) array_like
        The denominator coefficient vector of the filter.  If ``a[0]``
        is not 1, then both `a` and `b` are normalized by ``a[0]``.
    x : array_like
        The array of data to be filtered.
    axis : int, optional
        The axis of `x` to which the filter is applied.
        Default is -1.
    padtype : str or None, optional
        Must be 'odd', 'even', 'constant', or None.  This determines the
        type of extension to use for the padded signal to which the filter
        is applied.  If `padtype` is None, no padding is used.  The default
        is 'odd'.
    padlen : int or None, optional
        The number of elements by which to extend `x` at both ends of
        `axis` before applying the filter.  This value must be less than
        ``x.shape[axis] - 1``.  ``padlen=0`` implies no padding.
        The default value is ``3 * max(len(a), len(b))``.
    method : str, optional
        Determines the method for handling the edges of the signal, either
        "pad" or "gust".  When `method` is "pad", the signal is padded; the
        type of padding is determined by `padtype` and `padlen`, and `irlen`
        is ignored.  When `method` is "gust", Gustafsson's method is used,
        and `padtype` and `padlen` are ignored.
    irlen : int or None, optional
        When `method` is "gust", `irlen` specifies the length of the
        impulse response of the filter.  If `irlen` is None, no part
        of the impulse response is ignored.  For a long signal, specifying
        `irlen` can significantly improve the performance of the filter.

    Returns
    -------
    y : ndarray
        The filtered output with the same shape as `x`.

    See Also
    --------
    sosfiltfilt, lfilter_zi, lfilter, lfiltic, savgol_filter, sosfilt

    Notes
    -----
    When `method` is "pad", the function pads the data along the given axis
    in one of three ways: odd, even or constant.  The odd and even extensions
    have the corresponding symmetry about the end point of the data.  The
    constant extension extends the data with the values at the end points. On
    both the forward and backward passes, the initial condition of the
    filter is found by using `lfilter_zi` and scaling it by the end point of
    the extended data.

    When `method` is "gust", Gustafsson's method [1]_ is used.  Initial
    conditions are chosen for the forward and backward passes so that the
    forward-backward filter gives the same result as the backward-forward
    filter.

    The option to use Gustaffson's method was added in scipy version 0.16.0.

    References
    ----------
    .. [1] F. Gustaffson, "Determining the initial states in forward-backward
           filtering", Transactions on Signal Processing, Vol. 46, pp. 988-992,
           1996.

    Examples
    --------
    The examples will use several functions from `scipy.signal`.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    First we create a one second signal that is the sum of two pure sine
    waves, with frequencies 5 Hz and 250 Hz, sampled at 2000 Hz.

    >>> t = np.linspace(0, 1.0, 2001)
    >>> xlow = np.sin(2 * np.pi * 5 * t)
    >>> xhigh = np.sin(2 * np.pi * 250 * t)
    >>> x = xlow + xhigh

    Now create a lowpass Butterworth filter with a cutoff of 0.125 times
    the Nyquist frequency, or 125 Hz, and apply it to ``x`` with `filtfilt`.
    The result should be approximately ``xlow``, with no phase shift.

    >>> b, a = signal.butter(8, 0.125)
    >>> y = signal.filtfilt(b, a, x, padlen=150)
    >>> np.abs(y - xlow).max()
    9.1086182074789912e-06

    We get a fairly clean result for this artificial example because
    the odd extension is exact, and with the moderately long padding,
    the filter's transients have dissipated by the time the actual data
    is reached.  In general, transient effects at the edges are
    unavoidable.

    The following example demonstrates the option ``method="gust"``.

    First, create a filter.

    >>> b, a = signal.ellip(4, 0.01, 120, 0.125)  # Filter to be applied.
    >>> np.random.seed(123456)

    `sig` is a random input signal to be filtered.

    >>> n = 60
    >>> sig = np.random.randn(n)**3 + 3*np.random.randn(n).cumsum()

    Apply `filtfilt` to `sig`, once using the Gustafsson method, and
    once using padding, and plot the results for comparison.

    >>> fgust = signal.filtfilt(b, a, sig, method="gust")
    >>> fpad = signal.filtfilt(b, a, sig, padlen=50)
    >>> plt.plot(sig, 'k-', label='input')
    >>> plt.plot(fgust, 'b-', linewidth=4, label='gust')
    >>> plt.plot(fpad, 'c-', linewidth=1.5, label='pad')
    >>> plt.legend(loc='best')
    >>> plt.show()

    The `irlen` argument can be used to improve the performance
    of Gustafsson's method.

    Estimate the impulse response length of the filter.

    >>> z, p, k = signal.tf2zpk(b, a)
    >>> eps = 1e-9
    >>> r = np.max(np.abs(p))
    >>> approx_impulse_len = int(np.ceil(np.log(eps) / np.log(r)))
    >>> approx_impulse_len
    137

    Apply the filter to a longer signal, with and without the `irlen`
    argument.  The difference between `y1` and `y2` is small.  For long
    signals, using `irlen` gives a significant performance improvement.

    >>> x = np.random.randn(5000)
    >>> y1 = signal.filtfilt(b, a, x, method='gust')
    >>> y2 = signal.filtfilt(b, a, x, method='gust', irlen=approx_impulse_len)
    >>> print(np.max(np.abs(y1 - y2)))
    1.80056858312e-10

    """
    b = np.atleast_1d(b)
    a = np.atleast_1d(a)
    x = np.asarray(x)*/

    pair<int, Vectord> p = _validate_pad(padtype, padlen, x, max(a.size(), b.size()));
    int edge = p.first;
    Vectord ext = p.second;

    // # Get the steady state of the filter's step response.
    auto zi = lfilter_zi(b, a);

//    # Reshape zi and create x0 so that zi*x0 broadcasts
//    # to the correct value for the 'zi' keyword argument
//    # to lfilter.
//    zi_shape = [1] * x.ndim
//    zi_shape[axis] = zi.size
//    zi = np.reshape(zi, zi_shape)
    double x0 = ext[0]; // axis_slice(ext, stop=1, axis=axis)

//    # Forward filter.
    Vectord y = lfilter(b, a, ext, zi * x0);

//    # Backward filter.
//    # Create y0 so zi*y0 broadcasts appropriately.
    double y0 = y(Eigen::last);
    y = lfilter(b, a, y.reverse(), zi * y0);

//    # Reverse y.
    y.reverseInPlace();

    if (edge > 0)
//        # Slice the actual signal from the extended signal.
        y = y(Eigen::seqN(edge, x.size())); // stop=-edge, axis=axis)

    return y;
}
