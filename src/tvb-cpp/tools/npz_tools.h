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

#ifndef TVB_CPP_NPZ_TOOLS_H
#define TVB_CPP_NPZ_TOOLS_H

#include <functional>

#include <tvb-cpp/definitions.h>
#include <tvb-cpp/tools/cnpy.h>

namespace tvb {

    // Load from npz

    TArray2d_uptr npz2Matrixd(const std::string &filename, const std::string &index);
    TArray2d_uptr npy2Matrixd(const std::string &filename);

    TArray2dMap npz2MatrixdMap(const std::string &filename);

    std::vector<TArray2d> np2VecMatrixd(const std::string &filename, const std::string &index = "");

    std::vector<double>  npz2VecDouble(const std::string &filename, const std::string &index);

    // Save to npz

    void vecMatrixd2npz(const std::vector<TArray2d>& data, const std::string& filename, const std::string& index);
    void TArray1d2npz(const TArray1d& data, const std::string& filename, const std::string& index);

    void vecDouble2npz(const std::vector<double>& data, const std::string& filename, const std::string& index);

    void MatrixdMap2npz(const std::string &filename, const TArray2dMap& mmap);

    void Matrixd2np(const TArray2d& data, const std::string& filename, const std::string& index = "");


}

#endif //TVB_CPP_NPZ_TOOLS_H
