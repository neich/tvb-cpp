//
// Created by imartin on 28-Oct-20.
//

#ifndef TVB_CPP_SIGNALTOOLS_H
#define TVB_CPP_SIGNALTOOLS_H

#include <definitions.h>

tvb::Vectord detrend_linear(const tvb::Vectord& data);

tvb::Vectord filtfilt_pad(const tvb::Vectord& b, const tvb::Vectord& a, const tvb::Vectord& x,
                     int padlen = 0, const std::string& padtype = "odd");

#endif //TVB_CPP_SIGNALTOOLS_H
