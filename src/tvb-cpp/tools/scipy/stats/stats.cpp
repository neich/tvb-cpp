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

#include <unordered_map>
#include <vector>
#include <algorithm>

#include "tvb-cpp/definitions.h"
#include "stats.h"
#include "ksstats.h"
#include <tvb-cpp/tools/numpy/numpy.h>

using namespace std;

int binomialCoefficients(int n, int k) {
    tvb::TArray1di C = tvb::TArray1di::Zero(k + 1);
    C[0] = 1;
    for (int i = 1; i <= n; i++) {
        for (int j = min(i, k); j > 0; j--)
            C[j] = C[j] + C[j - 1];
    }
    return C[k];
}

double corrcoef(const TArray1d &x, const TArray1d &y) {
    assert(x.size() == y.size());
    double xy_s = 0.0;
    double x_s = 0.0;
    double y_s = 0.0;
    double x2_s = 0.0;
    double y2_s = 0.0;
    unsigned n = x.size();
    for (unsigned i = 0; i < n; ++i) {
        double xi = x[i];
        double yi = y[i];
        xy_s += xi*yi;
        x_s += xi;
        y_s += yi;
        x2_s += xi*xi;
        y2_s += yi*yi;
    }
    double num = xy_s - (x_s * y_s/(double)n);
    double d = (x2_s - x_s*x_s/(double)n) * (y2_s - y_s*y_s/(double)n);
    if (d <= 0.0) return 0.0;
    double den = pow(d, 0.5);
    return num / den;
}

int gcd(int a, int b)
{
    // Everything divides 0
    if (a == 0)
        return b;
    if (b == 0)
        return a;

    // base case
    if (a == b)
        return a;

    // a is greater
    if (a > b)
        return gcd(a-b, b);
    return gcd(a, b-a);
}

double compute_prob_outside_square(int n, int h) {
//    """
//    Compute the proportion of paths that pass outside the two diagonal lines.
//    Parameters
//    ----------
//    n : integer
//        n > 0
//    h : integer
//        0 <= h <= n
//    Returns
//    -------
//    p : float
//        The proportion of paths that pass outside the lines x-y = +/-h.
//    """
//    # Compute Pr(D_{n,n} >= h/n)
//    # Prob = 2 * ( binom(2n, n-h) - binom(2n, n-2a) + binom(2n, n-3a) - ... )  / binom(2n, n)
//    # This formulation exhibits subtractive cancellation.
//    # Instead divide each term by binom(2n, n), then factor common terms
//    # and use a Horner-like algorithm
//    # P = 2 * A0 * (1 - A1*(1 - A2*(1 - A3*(1 - A4*(...)))))

    double P = 0.0;
    int k = int(floor(n / h));
    while (k >= 0) {
        double p1 = 1.0;
        // # Each of the Ai terms has numerator and denominator with h simple terms.
        for(int j = 0; j < h; ++j)
            p1 = (n - k * h - j) * p1 / (n + k * h + j + 1);
        P = p1 * (1.0 - P);
        k -= 1;
    }
    return 2 * P;
}

int count_paths_outside_method(int m, int n, int g, int h) {
//    """
//    Count the number of paths that pass outside the specified diagonal.
//    Parameters
//    ----------
//    m : integer
//        m > 0
//    n : integer
//        n > 0
//    g : integer
//        g is greatest common divisor of m and n
//    h : integer
//        0 <= h <= lcm(m,n)
//    Returns
//    -------
//    p : float
//        The number of paths that go low.
//        The calculation may overflow - check for a finite answer.
//    Raises
//    ------
//    FloatingPointError: Raised if the intermediate computation goes outside
//    the range of a float.
//    Notes
//    -----
//    Count the integer lattice paths from (0, 0) to (m, n), which at some
//    point (x, y) along the path, satisfy:
//      m*y <= n*x - h*g
//    The paths make steps of size +1 in either positive x or positive y directions.
//    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.
//    Hodges, J.L. Jr.,
//    "The Significance Probability of the Smirnov Two-Sample Test,"
//    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.
//    """
//    # Compute #paths which stay lower than x/m-y/n = h/lcm(m,n)
//    # B(x, y) = #{paths from (0,0) to (x,y) without previously crossing the boundary}
//    #         = binom(x, y) - #{paths which already reached the boundary}
//    # Multiply by the number of path extensions going from (x, y) to (m, n)
//    # Sum.
//
//    # Probability is symmetrical in m, n.  Computation below assumes m >= n.
    if (m < n)
        std::swap(m, n);
    int mg = m / g;
    int ng = n / g;

//    # Not every x needs to be considered.
//    # xj holds the list of x values to be checked.
//    # Wherever n*x/m + ng*h crosses an integer
    int lxj = n + (mg-h)/mg;
    vector<int> xj;
    for (unsigned j = 0; j < lxj; ++j)
        xj.push_back((h + mg * j + ng-1)/ng);
//    # B is an array just holding a few values of B(x,y), the ones needed.
//    # B[j] == B(x_j, j)
    if (lxj == 0)
        return binomialCoefficients(m + n, n);
    TArray1d B(lxj);
    B.setZero();
    B[0] = 1;
//    # Compute the B(x, y) terms
//    # The binomial coefficient is an integer, but special.binom() may return a float.
//    # Round it to the nearest integer.
    for (unsigned j = 1; j < lxj; ++j) {
        double Bj = round(binomialCoefficients(xj[j] + j, j));
        if (!isfinite(Bj))
            throw "FloatingPointError";
        for (unsigned i = 0; i < j; ++i) {
            double bin = round(binomialCoefficients(xj[j] - xj[i] + j - i, j - i));
            Bj -= bin * B[i];
        }
        B[j] = Bj;
        if (!isfinite(Bj))
            throw "FloatingPointError";
    }
    // # Compute the number of path extensions...
    double num_paths = 0.0;
    for (unsigned j = 0; j < lxj; ++j) {
        double bin = round(binomialCoefficients((m - xj[j]) + (n - j), n - j));
        double term = B[j] * bin;
        if (!isfinite(term))
            throw "FloatingPointError";
        num_paths += term;
    }
    return int(round(num_paths));
}

double compute_prob_inside_method(int m, int n, int g, int h) {
//    """
//    Count the proportion of paths that stay strictly inside two diagonal lines.
//    Parameters
//    ----------
//    m : integer
//        m > 0
//    n : integer
//        n > 0
//    g : integer
//        g is greatest common divisor of m and n
//    h : integer
//        0 <= h <= lcm(m,n)
//    Returns
//    -------
//    p : float
//        The proportion of paths that stay inside the two lines.
//    Count the integer lattice paths from (0, 0) to (m, n) which satisfy
//    |x/m - y/n| < h / lcm(m, n).
//    The paths make steps of size +1 in either positive x or positive y directions.
//    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.
//    Hodges, J.L. Jr.,
//    "The Significance Probability of the Smirnov Two-Sample Test,"
//    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.
//    """
//    # Probability is symmetrical in m, n.  Computation below uses m >= n.
    if (m < n)
        std::swap(m, n);
    int mg = m / g;
    int ng = n / g;

//    # Count the integer lattice paths from (0, 0) to (m, n) which satisfy
//    # |nx/g - my/g| < h.
//    # Compute matrix A such that:
//    #  A(x, 0) = A(0, y) = 1
//    #  A(x, y) = A(x, y-1) + A(x-1, y), for x,y>=1, except that
//    #  A(x, y) = 0 if |x/m - y/n|>= h
//    # Probability is A(m, n)/binom(m+n, n)
//    # Optimizations exist for m==n, m==n*p.
//    # Only need to preserve a single column of A, and only a sliding window of it.
//    # minj keeps track of the slide.
    int minj = 0;
    int maxj = std::min(int(ceil(double(h) / double(mg))), n + 1);
    int curlen = maxj - minj;
//    # Make a vector long enough to hold maximum window needed.
    int lenA = std::min(2 * maxj + 2, n + 1);
//    # This is an integer calculation, but the entries are essentially
//    # binomial coefficients, hence grow quickly.
//    # Scaling after each column is computed avoids dividing by a
//    # large binomial coefficent at the end, but is not sufficient to avoid
//    # the large dyanamic range which appears during the calculation.
//    # Instead we rescale based on the magnitude of the right most term in
//    # the column and keep track of an exponent separately and apply
//    # it at the end of the calculation.  Similarly when multiplying by
//    # the binomial coefficint
    TArray1d A(lenA);
//    # Initialize the first column
    A(Eigen::seq(minj, maxj)) = 1.0;
    int expnt = 0;
    for (unsigned i = 0; i < m + 1; ++i) {
//        # Generate the next column.
//        # First calculate the sliding window
        int lastminj = minj;
        int lastlen = curlen;
        minj = max(int(floor((ng * i - h) / mg)) + 1, 0);
        minj = min(minj, n);
        maxj = min(int(ceil((ng * i + h) / mg)), n + 1);
        if (maxj <= minj)
            return 0;
//        # Now fill in the values
        A(Eigen::seq(0, maxj - minj)) = cumsum(A(Eigen::seq(minj - lastminj, maxj - lastminj - 1)));
        curlen = maxj - minj;
        if (lastlen > curlen)
//            # Set some carried-over elements to 0
            A(Eigen::seq(maxj - minj, maxj - minj + (lastlen - curlen))) = 0;
//        # Rescale if the right most value is over 2**900
        double val = A[maxj - minj - 1];
        int valexpt;
        frexp(val, &valexpt);
        if (valexpt > 900.0) {
//            # Scaling to bring down to about 2**800 appears
//            # sufficient for sizes under 10000.
            valexpt -= 800.0;
            std::transform(A.begin(), A.end(), A.begin(),
                           [valexpt](double m) { return m*pow(2.0, -valexpt); });
            expnt += valexpt;
        }
    }

    double val = A[maxj - minj - 1];
//    # Now divide by the binomial (m+n)!/m!/n!
    for (unsigned i = 0; i < n+1; ++i) {
        val = (val * i) / (m + i);
        int valexpt;
        frexp(val, &valexpt);
        if (valexpt < -128)
            val = val * pow(2.0, -valexpt);
        expnt += valexpt;
    }
//    # Finally scale if needed.
    return val * pow(2.0, expnt);
}


tuple<bool, double, double> attempt_exact_2kssamp(int n1, int n2, int g, double d,
                                            const std::string& alternative) {
//    """Attempts to compute the exact 2sample probability.
//    n1, n2 are the sample sizes
//    g is the gcd(n1, n2)
//    d is the computed max difference in ECDFs
//    Returns (success, d, probability)
//    """
    double qnan = std::numeric_limits<double>::quiet_NaN();
    int lcm = (n1 / g) * n2;
    int h = int(d * lcm);
    d = h * 1.0 / lcm;
    if (h == 0)
        return {true, d, 1.0};
    bool saw_fp_error = false;
    double prob = qnan;
    if (alternative =="two-sided") {
        if (n1 == n2)
            prob = compute_prob_outside_square(n1, h);
        else
            prob = 1 - compute_prob_inside_method(n1, n2, g, h);
    }
    else {
        if (n1 == n2) {
//                # prob = binom(2n, n-h) / binom(2n, n)
//                # Evaluating in that form incurs roundoff errors
//                # from special.binom. Instead calculate directly
            prob = 1.0;
            for (unsigned j = 0; j < h; ++j)
                prob *= (n1 - j) / (n1 + j + 1.0);
        } else {
            int num_paths = count_paths_outside_method(n1, n2, g, h);
            int bin = binomialCoefficients(n1 + n2, n1);
            // if (!isfinite(bin) || !isfinite(num_paths) || num_paths > bin)
            if (num_paths > bin)
                saw_fp_error = true;
            else
                prob = num_paths / bin;
        }
    }
    if (saw_fp_error)
        return { false, d, qnan };
    if (prob < 0.0 || prob > 1.0)
        return {false, d, prob};
    return {true, d, prob};
}

std::pair<double, double> ks_2samp(const TArray1d &d1,
                                   const TArray1d &d2,
                                   const std::string& alternative,
                                   std::string mode) {
    static vector<string> modes = {"auto", "exact", "asymp"};
    assert(("ks_2samp must take one of these three modes: auto, exact,asymp", find(modes.begin(), modes.end(), mode) != modes.end()));
    static vector<string> alternatives = {"two-sided", "greater", "less"};
    assert(("ks_2samp must take one of these three alternatives: two-sided, greater, less", find(alternatives.begin(), alternatives.end(), alternative) != alternatives.end()));
    unsigned MAX_AUTO_N = 10000;  // "auto" will attempt to be exact if n1,n2 <= MAX_AUTO_N
    TArray1d_sptr data1 = std::make_shared<TArray1d >(d1);
    TArray1d_sptr data2 = std::make_shared<TArray1d >(d2);
    sort(data1->begin(), data1->end());
    sort(data2->begin(), data2->end());
    auto n1 = data1->size();
    auto n2 = data2->size();
    assert(n1 > 0 && n2 > 0);

    TArray1d_sptr data_all = std::make_shared<TArray1d>(n1 + n2);
    *data_all << *data1, *data2;
    data1.reset();
    data2.reset();
    // using searchsorted solves equal data problem
    TArray1d_sptr cdf1 = std::make_shared<TArray1d >(searchsorted(data_all->block(0, 0, n1, 1), *data_all, "right").cast<Float>() / Float(n1));
    TArray1d_sptr cdf2 = std::make_shared<TArray1d >(searchsorted(data_all->block(n1, 0, n2, 1), *data_all, "right").cast<Float>() / Float(n2));
    data_all.reset();
//    TArray1d cdf1(cdfi1.size()), cdf2(cdfi2.size());
//    std::transform(cdfi1.begin(), cdfi1.end(), cdf1.begin(), [n1](int v) { return double(v)/n1; });
//    std::transform(cdfi2.begin(), cdfi2.end(), cdf2.begin(), [n2](int v) { return double(v)/n2; });

    TArray1d cddiffs = *cdf1 - *cdf2;
    cdf1.reset();
    cdf2.reset();
    Float minS = std::clamp(-cddiffs.minCoeff(), Float(0.0), Float(1.0));  // Ensure sign of minS is not negative.
    Float maxS = cddiffs.maxCoeff();
    unordered_map<string, Float>  alt2Dvalue = {{"less", minS},
                                                 {"greater", maxS},
                                                 {"two-sided", max(minS, maxS)}};
    Float d = alt2Dvalue[alternative];
    int g = gcd(n1, n2);
    int n1g = n1 / g;
    int n2g = n2 / g;
    double prob = -std::numeric_limits<double>::infinity();
    if (mode =="auto")
        mode = max(n1, n2) <= MAX_AUTO_N ? "exact" : "asymp";
    else if (mode =="exact") {
        // If lcm(n1, n2) is too big, switch from exact to asymp
        if (n1g >= std::numeric_limits<int>::max() / n2g)
            mode = "asymp";
//        warnings.warn(
//                f
//        "Exact ks_2samp calculation not possible with samples sizes "
//        f
//        "{n1} and {n2}. Switching to "asymp".", RuntimeWarning)
    }
    if (mode == "exact") {
        auto [success, dd, prob] = attempt_exact_2kssamp(n1, n2, g, d, alternative);
        d = dd;
        if (!success)
            mode = "asymp";
//        if original_mode == "exact":
//        warnings.warn(f
//        "ks_2samp: Exact calculation unsuccessful. "
//        f
//        "Switching to mode={mode}.", RuntimeWarning)
    }
    if (mode == "asymp") {
//# The product n1*n2 is large.  Use Smirnov"s asymptoptic formula.
//# Ensure float to avoid overflow in multiplication
//# sorted because the one-sided formula is not symmetric in n1, n2
        double m = n1;
        double n = n2;
        if (m < n) std::swap(m, n);
        double en = m * n / (m + n);
        if (alternative == "two-sided")
            prob = kolmogn(round(en), d, false);
        else {
            double z = sqrt(en) * d;
//# Use Hodges" suggested approximation Eqn 5.3
//# Requires m to be the larger of (n1, n2)
            double expt = -2 * z * z - 2 * z * (m + 2 * n) / sqrt(m * n * (m + n)) / 3.0;
            prob = exp(expt);
        }
    }
    prob = std::clamp(prob, 0.0, 1.0);
    return {d, prob};
}

double ks2_sf(double n, double x) {
    static double loc = 0.0;
    static double scale = 1.0;

    double _a = 0.5 / n;
    double _b = 1.0;
    x = (x - loc) /  scale;
    bool cond0 = n > 0.0;
    bool cond1 = (_a < x) && (x < _b);
    bool cond2 = cond0 && (x <= _a);
    bool cond = cond0 && cond1;
    double output = cond2 ? 1.0 : 0.0;
    return output;
}