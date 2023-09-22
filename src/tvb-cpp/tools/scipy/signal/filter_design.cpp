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

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <unordered_map>

#include "tvb-cpp/definitions.h"

#include <tvb-cpp/tools/numpy/numpy.h>
#include <tvb-cpp/tools/numpy/core/numeric.h>
#include <tvb-cpp/tools/eigen/eigen.h>

#include "filter_design.h"

using namespace std;
using namespace tvb;


static unordered_map<string, string> band_dict = {
        {"band", "bandpass"},
        {"bandpass", "bandpass"},
        {"pass", "bandpass"},
        {"bp", "bandpass"},

        {"bs", "bandstop"},
        {"bandstop", "bandstop"},
        {"bands", "bandstop"},
        {"stop", "bandstop"},

        {"l", "lowpass"},
        {"low", "lowpass"},
        {"lowpass", "lowpass"},
        {"lp", "lowpass"},

        {"high", "highpass"},
        {"highpass", "highpass"},
        {"h", "highpass"},
        {"hp", "highpass"},
};

int _relative_degree(const TArray1d& z, const TArray1dc& p) {
//    """
//    Return relative degree of transfer function from zeros and poles
//    """
    int degree = p.size() - z.size();
    if (degree < 0)
        throw("Improper transfer function. Must have at least as many poles as zeros.");
    else
        return degree;
}

ZeroPoleGain buttap(double N) {
//    """Return (z,p,k) for analog prototype of Nth-order Butterworth filter.
//
//    The filter will have an angular (e.g., rad/s) cutoff frequency of 1.
//
//    See Also
//    --------
//    butter : Filter design function using this prototype
//
//    """
    if (abs(int(N)) != N)
        throw("Filter order must be a nonnegative integer");
    TArray1d z;
    TArray1d m = arange<Float>(-N + 1, N, 2);
    int size = m.size();
    TArray1dc p(size);
    for (unsigned i = 0; i < size; ++i)
        p[i] = -exp(std::complex<tvb::Float>(1i) * tvb::Float(M_PI) * m[i] / (tvb::Float(2.0)*tvb::Float(size)));
    // # Middle value is 0 to ensure an exactly real pole
    // p = -numpy.exp(1i * pi * m / (2 * N))
    double k = 1;
    return {z, p, k};
}


ZeroPoleGain lp2lp_zpk(const TArray1d& z, const TArray1dc& p, double k, double wo=1.0) {
/*
    r"""
    Transform a lowpass filter prototype to a different frequency.

    Return an analog low-pass filter with cutoff frequency `wo`
    from an analog low-pass filter prototype with unity cutoff frequency,
    using zeros, poles, and gain ('zpk') representation.

    Parameters
    ----------
    z : array_like
        Zeros of the analog filter transfer function.
    p : array_like
        Poles of the analog filter transfer function.
    k : float
        System gain of the analog filter transfer function.
    wo : float
        Desired cutoff, as angular frequency (e.g., rad/s).
        Defaults to no change.

    Returns
    -------
    z : ndarray
        Zeros of the transformed low-pass filter transfer function.
    p : ndarray
        Poles of the transformed low-pass filter transfer function.
    k : float
        System gain of the transformed low-pass filter.

    See Also
    --------
    lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
    lp2lp

    Notes
    -----
    This is derived from the s-plane substitution

    .. math:: s \rightarrow \frac{s}{\omega_0}

    .. versionadded:: 1.1.0

    """
*/
//    z = atleast_1d(z)
//    p = atleast_1d(p)
    wo = float(wo); //  # Avoid int wraparound

    int degree = _relative_degree(z, p);

    // # Scale all points radially from origin to shift cutoff frequency
    TArray1d z_lp = wo * z;
    TArray1dc p_lp = wo * p;

//    # Each shifted pole decreases gain by wo, each shifted zero increases it.
//    # Cancel out the net change to keep overall gain the same
    double k_lp = k * pow(wo, degree);

    return {z_lp, p_lp, k_lp};
}

ZeroPoleGain lp2hp_zpk(const TArray1d& z, const TArray1dc& p, double k, double wo=1.0) {
//    r"""
//    Transform a lowpass filter prototype to a highpass filter.
//
//    Return an analog high-pass filter with cutoff frequency `wo`
//    from an analog low-pass filter prototype with unity cutoff frequency,
//    using zeros, poles, and gain ('zpk') representation.
//
//    Parameters
//    ----------
//    z : array_like
//        Zeros of the analog filter transfer function.
//    p : array_like
//        Poles of the analog filter transfer function.
//    k : float
//        System gain of the analog filter transfer function.
//    wo : float
//        Desired cutoff, as angular frequency (e.g., rad/s).
//        Defaults to no change.
//
//    Returns
//    -------
//    z : ndarray
//        Zeros of the transformed high-pass filter transfer function.
//    p : ndarray
//        Poles of the transformed high-pass filter transfer function.
//    k : float
//        System gain of the transformed high-pass filter.
//
//    See Also
//    --------
//    lp2lp_zpk, lp2bp_zpk, lp2bs_zpk, bilinear
//    lp2hp
//
//    Notes
//    -----
//    This is derived from the s-plane substitution
//
//    .. math:: s \rightarrow \frac{\omega_0}{s}
//
//    This maintains symmetry of the lowpass and highpass responses on a
//    logarithmic scale.
//
//    .. versionadded:: 1.1.0
//
//    """
//    z = atleast_1d(z)
//    p = atleast_1d(p)

    int degree = _relative_degree(z, p);

//    # Invert positions radially about unit circle to convert LPF to HPF
//    # Scale all points radially from origin to shift cutoff frequency
    TArray1d zz_hp = wo / z;
    TArray1dc p_hp = wo / p;

    // # If lowpass had zeros at infinity, inverting moves them to origin.
    TArray1d z_hp(zz_hp.size() + degree);
    z_hp << zz_hp, TArray1d::Zero(degree);

    // # Cancel out gain change caused by inversion
    double k_hp = k * real((-z).prod() / (-p).prod());

    return {z_hp, p_hp, k_hp};
}


ZeroPoleGain lp2bp_zpk(const TArray1d& z, const TArray1dc& p, double k,
                       double wo=1.0, double bw=1.0) {
/*
    r"""
    Transform a lowpass filter prototype to a bandpass filter.

    Return an analog band-pass filter with center frequency `wo` and
    bandwidth `bw` from an analog low-pass filter prototype with unity
    cutoff frequency, using zeros, poles, and gain ('zpk') representation.

    Parameters
    ----------
    z : array_like
        Zeros of the analog filter transfer function.
    p : array_like
        Poles of the analog filter transfer function.
    k : float
        System gain of the analog filter transfer function.
    wo : float
        Desired passband center, as angular frequency (e.g., rad/s).
        Defaults to no change.
    bw : float
        Desired passband width, as angular frequency (e.g., rad/s).
        Defaults to 1.

    Returns
    -------
    z : ndarray
        Zeros of the transformed band-pass filter transfer function.
    p : ndarray
        Poles of the transformed band-pass filter transfer function.
    k : float
        System gain of the transformed band-pass filter.

    See Also
    --------
    lp2lp_zpk, lp2hp_zpk, lp2bs_zpk, bilinear
    lp2bp

    Notes
    -----
    This is derived from the s-plane substitution

    .. math:: s \rightarrow \frac{s^2 + {\omega_0}^2}{s \cdot \mathrm{BW}}

    This is the "wideband" transformation, producing a passband with
    geometric (log frequency) symmetry about `wo`.

    .. versionadded:: 1.1.0

    """
*/
//    z = atleast_1d(z)
//    p = atleast_1d(p)
//    wo = float(wo)
//    bw = float(bw)

    int degree = _relative_degree(z, p);

    // # Scale poles and zeros to desired bandwidth
    TArray1d zz_lp = z * bw / 2.0;
    TArray1dc z_lp = zz_lp;
    TArray1dc p_lp = p * bw / 2.0;

    // # Square root needs to produce complex result, not NaN
//    z_lp = z_lp.astype(complex);
//    p_lp = p_lp.astype(complex);

    // # Duplicate poles and zeros and shift from baseband to +wo and -wo
    TArray1dc zz_bp(z_lp.size() * 2);
    TArray1dc p_bp(p_lp.size() * 2);
    zz_bp << z_lp + (z_lp.pow(2.0) - pow(wo, 2.0)).sqrt(),
            z_lp - (z_lp.pow(2.0) - pow(wo, 2.0)).sqrt();
    p_bp << p_lp + (p_lp.pow(2.0) - pow(wo, 2.0)).sqrt(),
            p_lp - (p_lp.pow(2.0) - pow(wo, 2.0)).sqrt();

    // # Move degree zeros to origin, leaving degree zeros at infinity for BPF
    TArray1dc z_bp(zz_bp.size() + degree);
    z_bp << zz_bp, TArray1dc::Zero(degree);

    // # Cancel out gain change from frequency scaling
    double k_bp = k * pow(bw, degree);

    return {vc2vd(z_bp), p_bp, k_bp};
}


ZeroPoleGain lp2bs_zpk(const TArray1d& z, const TArray1dc& p, double k,
                       Float wo=1.0, Float bw=1.0) {
//    r"""
//    Transform a lowpass filter prototype to a bandstop filter.
//
//    Return an analog band-stop filter with center frequency `wo` and
//    stopband width `bw` from an analog low-pass filter prototype with unity
//    cutoff frequency, using zeros, poles, and gain ('zpk') representation.
//
//    Parameters
//    ----------
//    z : array_like
//        Zeros of the analog filter transfer function.
//    p : array_like
//        Poles of the analog filter transfer function.
//    k : float
//        System gain of the analog filter transfer function.
//    wo : float
//        Desired stopband center, as angular frequency (e.g., rad/s).
//        Defaults to no change.
//    bw : float
//        Desired stopband width, as angular frequency (e.g., rad/s).
//        Defaults to 1.
//
//    Returns
//    -------
//    z : ndarray
//        Zeros of the transformed band-stop filter transfer function.
//    p : ndarray
//        Poles of the transformed band-stop filter transfer function.
//    k : float
//        System gain of the transformed band-stop filter.
//
//    See Also
//    --------
//    lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, bilinear
//    lp2bs
//
//    Notes
//    -----
//    This is derived from the s-plane substitution
//
//    .. math:: s \rightarrow \frac{s \cdot \mathrm{BW}}{s^2 + {\omega_0}^2}
//
//    This is the "wideband" transformation, producing a stopband with
//    geometric (log frequency) symmetry about `wo`.
//
//    .. versionadded:: 1.1.0
//
//    """
//    z = atleast_1d(z)
//    p = atleast_1d(p)
//    wo = float(wo)
//    bw = float(bw)

    int degree = _relative_degree(z, p);

    // # Invert to a highpass filter with desired bandwidth
    TArray1d z_hp = (bw / 2) / z;
    TArray1dc p_hp = (bw / 2) / p;

//    # Square root needs to produce complex result, not NaN
//    z_hp = z_hp.astype(complex)
//    p_hp = p_hp.astype(complex)
    TArray1d zz_lp = z * bw / 2.0;
    TArray1dc z_lp = zz_lp;
    TArray1dc p_lp = p * bw / 2.0;

    // # Duplicate poles and zeros and shift from baseband to +wo and -wo
    TArray1dc zz_bp(z_lp.size() * 2);
    TArray1dc p_bp(p_lp.size() * 2);
    zz_bp << z_lp + (z_lp.pow(2.0) - pow(wo, 2.0)).sqrt(),
            z_lp - (z_lp.pow(2.0) - pow(wo, 2.0)).sqrt();
    p_bp << p_lp + (p_lp.pow(2.0) - pow(wo, 2.0)).sqrt(),
            p_lp - (p_lp.pow(2.0) - pow(wo, 2.0)).sqrt();

    // # Move any zeros that were at infinity to the center of the stopband
    TArray1dc z_bp(zz_bp.size() + 2 * degree);
    z_bp << z_bp, TArray1dc::Constant(degree, std::complex<Float>(+1i) * wo), TArray1dc::Constant(degree, std::complex<Float>(-1i) * wo);

    // # Cancel out gain change caused by inversion
    double k_bp = k * pow(bw, degree);

    return {vc2vd(z_bp), p_bp, k_bp};
}

ZeroPoleGain bilinear_zpk(const TArray1d& z, const TArray1dc& p, double k, double fs) {
//    r"""
//    Return a digital IIR filter from an analog one using a bilinear transform.
//
//    Transform a set of poles and zeros from the analog s-plane to the digital
//    z-plane using Tustin's method, which substitutes ``(z-1) / (z+1)`` for
//    ``s``, maintaining the shape of the frequency response.
//
//    Parameters
//    ----------
//    z : array_like
//        Zeros of the analog filter transfer function.
//    p : array_like
//        Poles of the analog filter transfer function.
//    k : float
//        System gain of the analog filter transfer function.
//    fs : float
//        Sample rate, as ordinary frequency (e.g., hertz). No prewarping is
//        done in this function.
//
//    Returns
//    -------
//    z : ndarray
//        Zeros of the transformed digital filter transfer function.
//    p : ndarray
//        Poles of the transformed digital filter transfer function.
//    k : float
//        System gain of the transformed digital filter.
//
//    See Also
//    --------
//    lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk
//    bilinear
//
//    Notes
//    -----
//    .. versionadded:: 1.1.0
//
//    Examples
//    --------
//    >>> from scipy import signal
//    >>> import matplotlib.pyplot as plt
//
//    >>> fs = 100
//    >>> bf = 2 * np.pi * np.array([7, 13])
//    >>> filts = signal.lti(*signal.butter(4, bf, btype='bandpass', analog=True, output='zpk'))
//    >>> filtz = signal.lti(*signal.bilinear_zpk(filts.zeros, filts.poles, filts.gain, fs))
//    >>> wz, hz = signal.freqz_zpk(filtz.zeros, filtz.poles, filtz.gain)
//    >>> ws, hs = signal.freqs_zpk(filts.zeros, filts.poles, filts.gain, worN=fs*wz)
//    >>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hz).clip(1e-15)), label=r'$|H(j \omega)|$')
//    >>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hs).clip(1e-15)), label=r'$|H_z(e^{j \omega})|$')
//    >>> plt.legend()
//    >>> plt.xlabel('Frequency [Hz]')
//    >>> plt.ylabel('Magnitude [dB]')
//    >>> plt.grid()
//    """
//    z = atleast_1d(z)
//    p = atleast_1d(p)

    int degree = _relative_degree(z, p);

    double fs2 = 2.0 * fs;

    // # Bilinear transform the poles and zeros
    TArray1d z_z = TArray1d(fs2 + z) / TArray1d(fs2 - z);
    TArray1dc p_z = (fs2 + p) / (fs2 - p);

    // # Any zeros that were at infinity get moved to the Nyquist frequency
    TArray1d zz(z_z.size() + degree);
    zz << z_z, -TArray1d::Ones(degree);

    // # Compensate for gain change
    double k_z = k * real((fs2 - z).prod() / (fs2 - p).prod());

    return {zz, p_z, k_z};
}


pair<TArray1d, TArray1d> zpk2tf(const TArray1d& z, const TArray1dc& p, double k) {
/*
    """
    Return polynomial transfer function representation from zeros and poles

    Parameters
    ----------
    z : array_like
        Zeros of the transfer function.
    p : array_like
        Poles of the transfer function.
    k : float
        System gain.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    """
*/
//    z = atleast_1d(z)
//    k = atleast_1d(k)
//    if len(z.shape) > 1:
//        temp = poly(z[0])
//        b = zeros((z.shape[0], z.shape[1] + 1), temp.dtype.char)
//        if len(k) == 1:
//            k = [k[0]] * z.shape[0]
//        for i in range(z.shape[0]):
//            b[i] = k[i] * poly(z[i])
//    else:
    TArray1d b = k * poly(z);
    TArray1d a = poly(p);

//    # Use real output if possible. Copied from numpy.poly, since
//    # we can't depend on a specific version of numpy.
//    if issubclass(b.dtype.type, numpy.complexfloating):
//        //# if complex roots are all complex conjugates, the roots are real.
//        roots = numpy.asarray(z, complex)
//        pos_roots = numpy.compress(roots.imag > 0, roots)
//        neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
//        if len(pos_roots) == len(neg_roots):
//            if numpy.all(numpy.sort_complex(neg_roots) ==
//                         numpy.sort_complex(pos_roots)):
//                b = b.real.copy()
//
//    if issubclass(a.dtype.type, numpy.complexfloating):
//        //# if complex roots are all complex conjugates, the roots are real.
//        roots = numpy.asarray(p, complex)
//        pos_roots = numpy.compress(roots.imag > 0, roots)
//        neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
//        if len(pos_roots) == len(neg_roots):
//            if numpy.all(numpy.sort_complex(neg_roots) ==
//                         numpy.sort_complex(pos_roots)):
//                a = a.real.copy()

    return {b, a};
}


ZeroPoleGain _iirfilter(int N,
                           const tvb::TArray1d& wn,
                           float rp,
                           float rs,
                           const std::string& btype,
                           bool analog,
                           const std::string& ftype,
                           float fs) {


    // ftype, btype, output = [x.lower() for x in (ftype, btype, output)]
    // Wn = asarray(Wn)
    TArray1d Wn = wn;
    if (fs != 0.0) {
//        if (analog)
//            raise ValueError("fs cannot be specified for an analog filter")
        Wn = 2.0 * Wn / fs;
    }

    auto bt = band_dict[btype];

//    try:
//        typefunc = filter_dict[ftype][0]
//    except KeyError:
//        raise ValueError(""%s" is not a valid basic IIR filter." % ftype)

//    if output not in ["ba", "zpk", "sos"]:
//        raise ValueError(""%s" is not a valid output form." % output)
//
//    if rp is not None and rp < 0:
//        raise ValueError("passband ripple (rp) must be positive")
//
//    if rs is not None and rs < 0:
//        raise ValueError("stopband attenuation (rs) must be positive")
//
//    # Get analog lowpass prototype
    TArray1d z;
    TArray1dc p;
    double k;

    if (ftype == "butter")
        std::tie(z, p, k) = buttap(N);
//    elif typefunc == besselap:
//        z, p, k = typefunc(N, norm=bessel_norms[ftype])
//    elif typefunc == cheb1ap:
//        if rp is None:
//            raise ValueError("passband ripple (rp) must be provided to "
//                             "design a Chebyshev I filter.")
//        z, p, k = typefunc(N, rp)
//    elif typefunc == cheb2ap:
//        if rs is None:
//            raise ValueError("stopband attenuation (rs) must be provided to "
//                             "design an Chebyshev II filter.")
//        z, p, k = typefunc(N, rs)
//    elif typefunc == ellipap:
//        if rs is None or rp is None:
//            raise ValueError("Both rp and rs must be provided to design an "
//                             "elliptic filter.")
//        z, p, k = typefunc(N, rp, rs)
    else
        throw std::runtime_error(string_format("\"%s\" not implemented in iirfilter.", ftype.c_str()));

    // # Pre-warp frequencies for digital filter design
    TArray1d warped(Wn.size());
    if (!analog) {
        if ((Wn <= 0).any() || (Wn >= 1).any()) {
            if (fs > 0.0)
                throw std::runtime_error(string_format("Digital filter critical frequencies must be 0 < Wn < fs/2 (fs=%f -> fs/2=%f)", fs, fs/2));
            throw std::runtime_error("Digital filter critical frequencies must be 0 < Wn < 1");
            }
        fs = 2.0;
        for (unsigned i = 0; i < warped.size(); ++i)
            warped[i] = 2 * fs * tan(M_PI * Wn[i] / fs);
    }
    else
        warped = Wn;

    // # transform to lowpass, bandpass, highpass, or bandstop
    TArray1d zz;
    TArray1dc pp;
    double kk;
    if (btype == "lowpass" || btype == "highpass") {
        if (Wn.size() != 1)
            throw std::runtime_error("Must specify a single critical frequency Wn for lowpass or highpass filter");

        if (btype == "lowpass")
            std::tie(zz, pp, kk) = lp2lp_zpk(z, p, k, warped[0]);
        else if (btype == "highpass")
            std::tie(zz, pp, kk) = lp2hp_zpk(z, p, k, warped[0]);
    }
    else if (btype == "bandpass" || btype == "bandstop") {
        double bw, wo;
        try {
            bw = warped[1] - warped[0];
            wo = sqrt(warped[0] * warped[1]);
        }
        catch (...) {
            throw std::runtime_error("Wn must specify start and stop frequencies for bandpass or bandstop filter");
        }
        if (btype == "bandpass")
            std::tie(zz, pp, kk) = lp2bp_zpk(z, p, k, wo, bw);
        else if (btype == "bandstop")
            std::tie(zz, pp, kk) = lp2bs_zpk(z, p, k, wo, bw);
    }
    else
        throw std::runtime_error(string_format("\"%s\" not implemented in iirfilter.", btype.c_str()));

    //# Find discrete equivalent if necessary
    if (!analog)
        std::tie(zz, pp, kk) = bilinear_zpk(zz, pp, kk, fs);

    return {zz, pp, kk};

//    // # Transform to proper out type (pole-zero, state-space, numer-denom)
//    if (output == "zpk")
//        return {zz, pp, kk};
//    else if (output == "ba")
//        return zpk2tf(zz, pp, kk);
//    else if (output == "sos")
//        return zpk2sos(zz, pp, kk);
}

ZeroPoleGain iirfilter_zpk(int N,
                           const tvb::TArray1d& wn,
                           float rp,
                           float rs,
                           const std::string& btype,
                           bool analog,
                           const std::string& ftype,
                           float fs) {
    auto [z, p, k] = _iirfilter(N, wn, rp, rs, btype, analog, ftype, fs);
    return {z, p, k};
}

pair<TArray1d, TArray1d> iirfilter_ba(int N,
                                      const tvb::TArray1d& wn,
                                      float rp,
                                      float rs,
                                      const std::string& btype,
                                      bool analog,
                                      const std::string& ftype,
                                      float fs) {
    auto[z, p, k] = _iirfilter(N, wn, rp, rs, btype, analog, ftype, fs);
    return zpk2tf(z, p, k);
}
