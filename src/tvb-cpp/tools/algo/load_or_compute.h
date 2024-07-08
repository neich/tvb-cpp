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

#ifndef TVB_CPP_LOAD_OR_COMPUTE_H
#define TVB_CPP_LOAD_OR_COMPUTE_H

#include <fstream>

#include <tvb-cpp/definitions.h>
#include <tvb-cpp/tools/npz_tools.h>

namespace tvb {

    TArray2dMap
    load_or_compute(const std::string &filename,
                    const std::function<TArray2dMap()> &func) {
        std::ifstream f(filename.c_str());
        if (f.good())
            return npz2MatrixdMap(filename);
        else {
            TArray2dMap result = func();
            MatrixdMap2npz(filename, result);
            return result;
        }
    }

//    TArray1d
//    load_or_compute_index_1d(const std::string &filename, const std::string &index,
//                          const std::function<TArray1d()> &func) {
//        std::ifstream f(filename.c_str());
//        if (f.good())
//            return npz2Matrixd(filename, index);
//        else {
//            TArray1d result = func();
//            MatrixdMap2npz(result, filename, index);
//            return result;
//        }
//    }
//
//    TArray2d
//    load_or_compute_index_2d(const std::string &filename, const std::string &index,
//                          const std::function<TArray2d()> &func) {
//        std::ifstream f(filename.c_str());
//        if (f.good())
//            return npz2Matrixd(filename, index);
//        else {
//            TArray2d result = func();
//            Matrixd2npz(result, filename, index);
//            return result;
//        }
//    }


}

#endif //TVB_CPP_LOAD_OR_COMPUTE_H
