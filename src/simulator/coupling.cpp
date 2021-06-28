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

#include <simulator/coupling.h>

using namespace tvb;

Matrixd Coupling::couple(int step, const History &history) const {

    Matrixd a(3,3);
    Vectord b(3);

    AVectord c = a.matrix()*b.matrix();

    const Matrixd &g_ij = history.weights();
    History::QResult q = history.query(step);
    std::vector<Matrixd> pre = this->pre(q.first, q.second);
    const std::vector<int> c_vars = history.c_vars();
    Matrixd sum(history.num_nodes(), c_vars.size());
    for (int node = 0; node < history.num_nodes(); ++node) {
        for (int nv = 0; nv < c_vars.size(); ++nv) {
            Vectord v = pre[node].row(nv);
            sum(node, nv) = (g_ij.col(node) * v).sum();
        }
    }
    // Matrixd sum = (prod).rowwise().sum(); // sum!!!
    // Matrixd sum = (prod).colwise().sum().transpose(); // sum!!!
    return this->post(sum);
}