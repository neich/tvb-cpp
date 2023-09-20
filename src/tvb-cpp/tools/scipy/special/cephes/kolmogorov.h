//
// Created by imartin on 13-Oct-20.
//

#ifndef TVB_CPP_KOLMOGOROV_H
#define TVB_CPP_KOLMOGOROV_H


/* Struct to hold the CDF, SF and PDF, which are computed simultaneously */
typedef struct ThreeProbs {
    double sf;
    double cdf;
    double pdf;
} ThreeProbs;

ThreeProbs smirnov(int n, double x);
//double kolmogn(double n, double x, bool cdf = true);


#endif //TVB_CPP_KOLMOGOROV_H
