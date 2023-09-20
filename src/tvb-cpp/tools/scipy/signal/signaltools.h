//
// Created by imartin on 28-Oct-20.
//

#ifndef TVB_CPP_SIGNALTOOLS_H
#define TVB_CPP_SIGNALTOOLS_H

#include "tvb-cpp/definitions.h"

tvb::TArray1d detrend_linear(const tvb::TArray1d& data);

tvb::TArray1d filtfilt_pad(const tvb::TArray1d& b, const tvb::TArray1d& a, const tvb::TArray1d& x,
                           int padlen = 0, const std::string& padtype = "odd");

tvb::TArray2dc hilbert_matrix(const Eigen::Ref<const Eigen::MatrixXd> &bfFrame);

tvb::TArray1dc hilbert_signal(const tvb::TArray1d &x, int N = 0);


#endif //TVB_CPP_SIGNALTOOLS_H
