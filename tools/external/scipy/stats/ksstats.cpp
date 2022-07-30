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

#include <external/scipy/stats/ksstats.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

#include <definitions.h>

#include <external/scipy/special/cephes/kolmogorov.h>
#include <external/numpy/numpy.h>
#include <external/numpy/core/numeric.h>


using namespace tvb;
using namespace Eigen;
using namespace std;

const int E128 = 128;
const Float EP128 = pow(2, E128);
const Float EM128 = pow(2, -E128);

const Float _SQRT2PI = sqrt(2 * M_PI);
const Float _LOG_2PI = log(2 * M_PI);
const Float _MIN_LOG = -708;
const Float _SQRT3 = sqrt(3);
const Float _PI_SQUARED = pow(M_PI, 2);
const Float _PI_FOUR = pow(M_PI, 4);
const Float _PI_SIX = pow(M_PI, 6);


const Float LOG_2PI = log(2.0 * M_PI);
const std::vector<Float> STIRLING_COEFFS = {-2.955065359477124183e-2, 6.4102564102564102564e-3,
                                             -1.9175269175269175269e-3, 8.4175084175084175084e-4,
                                             -5.952380952380952381e-4, 7.9365079365079365079e-4,
                                             -2.7777777777777777778e-3, 8.3333333333333333333e-2};

int ipow(int base, int exp) {
    int result = 1;
    for (;;) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}


inline
TArray1d clip_prob(const TArray1d &p) {
    // """clips a probability to range 0<=p<=1."""
    TArray1d clipped(p.size());
    std::transform(p.begin(), p.end(), clipped.begin(),
                   [](Float v) { return std::clamp(v, Float(0.0), Float(1.0)); });
    return clipped;
}

inline
TArray1d select_and_clip_prob(const TArray1d &cdfprob, const TArray1d &sfprob, bool cdf = true) {
    // """Selects either the CDF or SF, and then clips to range 0<=p<=1."""
    const TArray1d &p = cdf ? cdfprob : sfprob;
    return clip_prob(p);
}

inline
Float select_and_clip_prob(Float cdfprob, Float sfprob, bool cdf = true) {
    // """Selects either the CDF or SF, and then clips to range 0<=p<=1."""
    Float p = cdf ? cdfprob : sfprob;
    return std::clamp(p, Float(0.0), Float(1.0));
}

Float log_nfactorial_div_n_pow_n(int n) {
//    # Computes n! / n**n
//    #    = (n-1)! / n**(n-1)
//    # Uses Stirling's approximation, but removes n*log(n) up-front to
//    # avoid subtractive cancellation.
//    #    = log(n)/2 - n + log(sqrt(2pi)) + sum B_{2j}/(2j)/(2j-1)/n**(2j-1)
    Float rn = 1.0 / n;
    return log(n) / tvb::Float(2.0) - n + LOG_2PI / tvb::Float(2.0) + rn * polyval(STIRLING_COEFFS, rn / n);
}


Float kolmogn_DMTW(int n, Float d, bool cdf = true) {
//    """Computes the Kolmogorov CDF:  Pr(D_n <= d) using the MTW approach to
//    the Durbin matrix algorithm.
//
//    Durbin (1968); Marsaglia, Tsang, Wang (2003). [1], [3].
//    """
//    # Write d = (k-h)/n, where k is positive integer and 0 <= h < 1
//    # Generate initial matrix H of size m*m where m=(2k-1)
//    # Compute k-th row of (n!/n^n) * H^n, scaling intermediate results.
//    # Requires memory O(m^2) and computation O(m^2 log(n)).
//    # Most suitable for small m.

    if (d >= 1.0)
        return select_and_clip_prob(1.0, 0.0, cdf);
    Float nd = n * d;
    if (nd <= 0.5)
        return select_and_clip_prob(0.0, 1.0, cdf);
    int k = int(ceil(nd));
    int h = k - nd;
    int m = 2 * k - 1;

    TMatrix H(m, m); // H = np.zeros([m, m])

//    # Initialize: v is first column (and last row) of H
//    #  v[j] = (1-h^(j+1)/(j+1)!  (except for v[-1])
//    #  w[j] = 1/(j)!
//    # q = k-th row of H (actually i!/n^i*H^i)
    // intm = np.arange(1, m + 1)
    // v = 1.0 - h ** intm
    TArray1d v(m);
    for (unsigned i = 0; i < m; ++i)
        v[i] = 1.0 - ipow(h, i + 1);
    TArray1d w = TArray1d::Zero(m);
    Float fac = 1.0;
    for (unsigned j = 1; j <= m; ++j) {
        w[j - 1] = fac;
        fac /= j; // # This might underflow.  Isn't a problem.
        v[j - 1] *= fac;
    }
    Float tt = pow(std::max(2 * h - 1, 0), m) - 2 * ipow(h, m);
    v[m - 1] = (1.0 + tt) * fac;

    for (unsigned i = 1; i <= m; ++i)
        H(seq(i - 1, last), i) = w(seq(0, m - i - 1));
    // H[i - 1:, i] = w[:m - i + 1]
    H(all, 0) = v;
    H(last, all) = v.reverse();

    TMatrix Hpwr = TMatrix::Identity(m, m);
    int nn = n;
    Float expnt = 0; //  # Scaling of Hpwr
    Float Hexpnt = 0; //  # Scaling of H
    while (nn > 0) {
        if (nn % 2 == 1) {
            Hpwr = Hpwr * H; // np.matmul(Hpwr, H)
            expnt += Hexpnt;
        }
        H = H * H;
        Hexpnt *= 2;
        // # Scale as needed.
        if (abs(H(k - 1, k - 1)) > EP128) {
            H /= EP128;
            Hexpnt += E128;
        }
        nn = nn / 2;
    }

    Float p = Hpwr(k - 1, k - 1);

    // # Multiply by n!/n^n
    for (unsigned i = 1; i < n + 1; ++i) {
        p = i * p / n;
        if (abs(p) < EM128) {
            p *= EP128;
            expnt -= E128;
        }
    }

    // # unscale
    if (expnt != 0)
        p = p * pow(2.0, expnt);

    return select_and_clip_prob(p, 1.0 - p, cdf);

}

pair<int, int> pomeranz_compute_j1j2(int i, int n, int ll, int ceilf, int roundf) {
    // """Compute the endpoints of the interval for row i."""
    int j1, j2;
    if (i == 0) {
        j1 = -ll - ceilf - 1;
        j2 = ll + ceilf - 1;
    } else {
        //# i + 1 = 2*ip1div2 + ip1mod2
        int ip1div2 = (i + 1 / 2);
        int ip1mod2 = (i + 1) % 2;
        if (ip1mod2 == 0) { // # i is odd
            if (ip1div2 == n + 1) {
                j1 = n - ll - ceilf - 1;
                j2 = n + ll + ceilf - 1;
            } else {
                j1 = ip1div2 - 1 - ll - roundf - 1;
                j2 = ip1div2 + ll - 1 + ceilf - 1;
            }
        } else {
            j1 = ip1div2 - 1 - ll - 1;
            j2 = ip1div2 + ll + roundf - 1;
        }
    }
    return {std::max(j1 + 2, 0), std::min(j2, n)};
}


Float kolmogn_Pomeranz(int n, Float x, bool cdf = true) {
//    r"""Computes Pr(D_n <= d) using the Pomeranz recursion algorithm.
//
//    Pomeranz (1974) [2]
//    """
//
//    # V is n*(2n+2) matrix.
//    # Each row is convolution of the previous row and probabilities from a
//    #  Poisson distribution.
//    # Desired CDF probability is n! V[n-1, 2n+1]  (final entry in final row).
//    # Only two rows are needed at any given stage:
//    #  - Call them V0 and V1.
//    #  - Swap each iteration
//    # Only a few (contiguous) entries in each row can be non-zero.
//    #  - Keep track of start and end (j1 and j2 below)
//    #  - V0s and V1s track the start in the two rows
//    # Scale intermediate results as needed.
//    # Only a few different Poisson distributions can occur
    Float t = n * x;
    int ll = int(floor(t));
    Float f = 1.0 * (t - ll); //  # fractional part of t
    Float g = std::min(f, tvb::Float(1.0) - f);
    int ceilf = f > 0 ? 1 : 0; //(1 if f > 0 else 0)
    int roundf = f > 0.5 ? 1 : 0; //(1 if f > 0.5 else 0)
    int npwrs = 2 * (ll + 1); //    # Maximum number of powers needed in convolutions
    TArray1d gpower = TArray1d::Zero(npwrs); // # gpower = (g/n)^m/m!
    TArray1d twogpower = TArray1d::Zero(npwrs); // # twogpower = (2g/n)^m/m!
    TArray1d onem2gpower = TArray1d::Zero(npwrs); // # onem2gpower = ((1-2g)/n)^m/m!
    // # gpower etc are *almost* Poisson probs, just missing normalizing factor.

    gpower[0] = 1.0;
    twogpower[0] = 1.0;
    onem2gpower[0] = 1.0;
    int expnt = 0;
    Float g_over_n = g / n;
    Float two_g_over_n = 2 * g / n;
    Float one_minus_two_g_over_n = (1 - 2 * g) / n;
    for (unsigned m = 1; m < npwrs; ++m) {
        gpower[m] = gpower[m - 1] * g_over_n / m;
        twogpower[m] = twogpower[m - 1] * two_g_over_n / m;
        onem2gpower[m] = onem2gpower[m - 1] * one_minus_two_g_over_n / m;
    }

    TArray1d V0 = TArray1d::Zero(npwrs);
    TArray1d V1 = TArray1d::Zero(npwrs);
    V1[0] = 1.0; //  # first row
    int V0s = 0;
    int V1s = 0; // # start indices of the two rows

    auto[j1, j2] = pomeranz_compute_j1j2(0, n, ll, ceilf, roundf);
    for (unsigned i = 1; i < 2 * n + 2; ++i) { //in range(1, 2 * n + 2):
        // # Preserve j1, V1, V1s, V0s from last iteration
        int k1 = j1;
        V0.swap(V1);
        swap(V0s, V1s);
        V1.setZero();
        auto[j1, j2] = pomeranz_compute_j1j2(i, n, ll, ceilf, roundf);
        TArray1d pwrs;
        if (i == 1 || i == 2 * n + 1)
            pwrs = gpower;
        else
            pwrs = i % 2 == 1 ? twogpower : onem2gpower;
        int ln2 = j2 - k1 + 1;
        if (ln2 > 0) {
            TArray1d conv = convolve(TArray1d(V0(seq(k1 - V0s, k1 - V0s + ln2 - 1))), TArray1d(pwrs(seq(0, ln2 - 1))));
            int conv_start = j1 - k1; //  # First index to use from conv
            int conv_len = j2 - j1 + 1; //  # Number of entries to use from conv
            V1(seq(0, conv_len - 1)) = conv(seq(conv_start, conv_start + conv_len - 1));
            // # Scale to avoid underflow.
            if (0 < V1.maxCoeff() < EM128) {
                V1 *= EP128;
                expnt -= E128;
            }
            V1s = V0s + j1 - k1;
        }
    }

    // # multiply by n!
    Float ans = V1[n - V1s];
    for (unsigned m = 1; m < n + 1; ++m) {
        if (abs(ans) > EP128) {
            ans *= EM128;
            expnt += E128;
        }
        ans *= m;
    }

    // # Undo any intermediate scaling
    if (expnt != 0)
        ans = ans * pow(2.0, expnt);
    ans = select_and_clip_prob(ans, 1.0 - ans, cdf);
    return ans;
}

Float kolmogn_PelzGood(int n, Float x, bool cdf = true) {
//    """Computes the Pelz-Good approximation to Prob(Dn <= x) with 0<=x<=1.
//
//    start with Li-Chien, Korolyuk approximation:
//        Prob(Dn <= x) ~ K0(z) + K1(z)/sqrt(n) + K2(z)/n + K3(z)/n**1.5
//    where z = x*sqrt(n).
//    Transform each K_(z) using Jacobi theta functions into a form suitable
//    for small z.
//    Pelz-Good (1976). [6]
//    """
    if (x <= 0.0)
        return select_and_clip_prob(0.0, 1.0, cdf = cdf);
    if (x >= 1.0)
        return select_and_clip_prob(1.0, 0.0, cdf = cdf);

    Float z = sqrt(n) * x;
    Float zsquared = pow(z, 2);
    Float zthree = pow(z, 3);
    Float zfour = pow(z, 4);
    Float zsix = pow(z, 6);

    Float qlog = -_PI_SQUARED / 8 / zsquared;
    if (qlog < _MIN_LOG) //  # z ~ 0.041743441416853426
        return select_and_clip_prob(0.0, 1.0, cdf = cdf);

    Float q = exp(qlog);

    // # Coefficients of terms in the sums for K1, K2 and K3
    Float k1a = -zsquared;
    Float k1b = _PI_SQUARED / 4;

    Float k2a = 6 * zsix + 2 * zfour;
    Float k2b = (2 * zfour - 5 * zsquared) * _PI_SQUARED / 4;
    Float k2c = _PI_FOUR * (1 - 2 * zsquared) / 16;

    Float k3d = _PI_SIX * (5 - 30 * zsquared) / 64;
    Float k3c = _PI_FOUR * (-60 * zsquared + 212 * zfour) / 16;
    Float k3b = _PI_SQUARED * (135 * zfour - 96 * zsix) / 4;
    Float k3a = -30 * zsix - 90 * pow(z, 8);

    TArray1d K0to3 = TArray1d::Zero(4);
//    # Use a Horner scheme to evaluate sum c_i q^(i^2)
//    # Reduces to a sum over odd integers.
    int maxk = int(ceil(16 * z / M_PI));
    for (unsigned k = maxk; k > 0; --k) {
        int m = 2 * k - 1;
        Float msquared = pow(m, 2);
        Float mfour = pow(m, 4);
        Float msix = pow(m, 6);
        Float qpower = pow(q, 8 * k);
        TArray1d coeffs(4);
        coeffs << 1.0,
                k1a + k1b * msquared,
                k2a + k2b * msquared + k2c * mfour,
                k3a + k3b * msquared + k3c * mfour + k3d * msix;
        K0to3 *= qpower;
        K0to3 += coeffs;
    }
    K0to3 *= q;
    K0to3 *= _SQRT2PI;
    // # z**10 > 0 as z > 0.04
    TArray1d div(4);
    div << z, 6 * zfour, 72 * pow(z, 7), 6480 * pow(z, 10);
    K0to3 /= div;

//    # Now do the other sum over the other terms, all integers k
//    # K_2:  (pi^2 k^2) q^(k^2),
//    # K_3:  (3pi^2 k^2 z^2 - pi^4 k^4)*q^(k^2)
//    # Don't expect much subtractive cancellation so use direct calculation
    q = exp(-_PI_SQUARED / 2 / zsquared);
    TArray1d ks(maxk);
    for (unsigned i = 0; i < maxk; ++i)
        ks[i] = maxk - i;  // np.arange(maxk, 0, -1)
    auto ksquared = ks.unaryExpr([](Float e) { return pow(Float(2.0), e); });
    Float sqrt3z = _SQRT3 * z;
    auto kspi = ks * M_PI;
    auto qpwers = ksquared.unaryExpr([q](Float k) { return pow(q, k); });
    Float k2extra = (ksquared * qpwers).sum();
    k2extra *= _PI_SQUARED * _SQRT2PI / (-36 * zthree);
    K0to3[2] += k2extra;
    Float k3extra = ((sqrt3z + kspi) * (sqrt3z - kspi) * ksquared * qpwers).sum();
    k3extra *= _PI_SQUARED * _SQRT2PI / (216 * zsix);
    K0to3[3] += k3extra;
    TArray1d powers_of_n = TArray1d::Constant(K0to3.size(), n * 1.0).pow(arange<Float>(0, K0to3.size()) / 2.0);
    K0to3 /= powers_of_n;

    if (!cdf) {
        K0to3 *= -1;
        K0to3[0] += 1;
    }

    Float Ksum = K0to3.sum();
    return Ksum;
}

Float kolmogn(Float n, Float x, bool cdf) {
//    """Computes the CDF(or SF) for the two-sided Kolmogorov-Smirnov statistic.
//
//    x must be of type float, n of type integer.
//
//    Simard & L'Ecuyer (2011) [7].
//    """
    if (x >= 1.0)
        return select_and_clip_prob(1.0, 0.0, cdf = cdf);
    if (x <= 0.0)
        return select_and_clip_prob(0.0, 1.0, cdf = cdf);
    Float t = n * x;
    Float prob = 1.0;
    if (t <= 1.0) { // # Ruben-Gambino: 1/2n <= x <= 1/n
        if (t <= 0.5)
            return select_and_clip_prob(0.0, 1.0, cdf = cdf);
        if (n <= 140) {
            for (unsigned i = 1; i < n + 1; ++i)
                prob *= i * (1.0 / n) * (2 * t - 1);
        } else
            prob = exp(log_nfactorial_div_n_pow_n(n) + n * (2 * t - 1));
        return select_and_clip_prob(prob, 1.0 - prob, cdf = cdf);
    }

    if (t >= n - 1) {//# Ruben-Gambino
        prob = 2 * pow((1.0 - x), n);
        return select_and_clip_prob(1 - prob, prob, cdf = cdf);
    }
    if (x >= 0.5) { //# Exact: 2 * smirnov
        prob = 2 * smirnov(n, x).sf;
        return select_and_clip_prob(1.0 - prob, prob, cdf = cdf);
    }

    Float nxsquared = t * x;
    if (n <= 140) {
        if (nxsquared <= 0.754693) {
            prob = kolmogn_DMTW(n, x, cdf = true);
            return select_and_clip_prob(prob, 1.0 - prob, cdf = cdf);
        }
        if (nxsquared <= 4) {
            prob = kolmogn_Pomeranz(n, x, cdf = true);
            return select_and_clip_prob(prob, 1.0 - prob, cdf = cdf);
        }
        // # Now use Miller approximation of 2*smirnov
        prob = 2 * smirnov(n, x).sf;
        return select_and_clip_prob(1.0 - prob, prob, cdf = cdf);
    }
    // # Split CDF and SF as they have different cutoffs on nxsquared.
    if (!cdf) {
        if (nxsquared >= 370.0)
            return 0.0;
        if (nxsquared >= 2.2) {
            int nn = (int)n;
            prob = 2.0 * smirnov(nn, x).sf;
            return std::clamp(prob, Float(0.0), Float(1.0));
        }
    }
    //# Fall through and compute the SF as 1.0-CDF
    Float cdfprob;
    if (nxsquared >= 18.0)
        cdfprob = 1.0;
    else if (n <= 100000 && n * pow(x, 1.5) <= 1.4)
        cdfprob = kolmogn_DMTW(n, x, cdf = true);
    else
        cdfprob = kolmogn_PelzGood(n, x, cdf = true);
    return select_and_clip_prob(cdfprob, 1.0 - cdfprob, cdf = cdf);
}
