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

#ifndef TVB_CPP_NUMERIC_H
#define TVB_CPP_NUMERIC_H

#include <definitions.h>

tvb::Vectord poly(const tvb::Vectord& seq_of_zeros);
tvb::Vectord poly(const tvb::AVectorc& seq_of_zeros);

tvb::Vectord convolve(const tvb::Vectord &a, const tvb::Vectord &b);
tvb::AVectorc convolve(const tvb::AVectorc &a, const tvb::AVectorc &b);


#endif //TVB_CPP_NUMERIC_H
