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

#ifndef TVB_CPP_CIRCSHIFT_H
#define TVB_CPP_CIRCSHIFT_H


#include <definitions.h>

using namespace tvb;

inline Matrixd circshift(Matrixd data, int a, int b = 0) {

    unsigned row = data.rows();
    unsigned col = data.cols();
    Matrixd out(row, col);
    Matrixd y(row, col);
    if (a > 0) // move down a line
    {
        out.topRows(a) = data.bottomRows(a);
        out.bottomRows(row - a) = data.topRows(row - a);
    }
    if (a < 0) // move up a line
    {
        out.topRows(row + a) = data.bottomRows(row + a);
        out.bottomRows(abs(a)) = data.topRows(abs(a));
    }
    if (a == 0) {
        out.array() = data.array();
    }

    if (b > 0) // move the b column to the right
    {
        y.leftCols(b) = out.rightCols(b);
        y.rightCols(col - b) = out.leftCols(col - b);
    }
    if (b < 0) // Move column b to the left
    {
        y.leftCols(col + b) = out.rightCols(col + b);
        y.rightCols(abs(b)) = out.leftCols(abs(b));

    }
    if (b == 0) {
        y = out;
    }
    return y;
}

#endif //TVB_CPP_CIRCSHIFT_H
