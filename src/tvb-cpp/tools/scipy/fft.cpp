//
// Created by natxm on 20/09/2023.
//

#include "fft.h"

using namespace tvb;

TArray1dc fft_inv(const TArray1dc &signal) {
    int N = signal.size();
    int i = N, j = N;
    const tvb::Float pi = 3.141592653589793238462643383279502884;

    TArray1dc table1 = signal;
    TArray1dc table3 = TArray1dc::Zero(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
//complex number method:
            table3[i] = std::complex<tvb::Float>(table3[i].real() + table1[j].real() * cos(i * j * 2 * pi / (tvb::Float) N),
                                                 table3[i].imag() + table1[j].real() * sin(i * j * 2 * pi / (tvb::Float) N));
            table3[i] = std::complex<tvb::Float>(table3[i].real() + table1[j].imag() * sin(i * j * 2 * pi / (tvb::Float) N) * -1,
                                                 table3[i].imag() + table1[j].imag() * cos(i * j * 2 * pi / (tvb::Float) N));

        }
    }
    table1 = table3 / N;
//        for (int j = 0; j < N; j++) {
//            table1[j].real() = std::complex<tvb::Float>(table3[j].real() / N, table3[j].imag() / N);
//        }

    return table1;
}

