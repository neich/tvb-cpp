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

#ifndef TVB_CPP_NOISE_H
#define TVB_CPP_NOISE_H

#include <random>

#include <definitions.h>

class UniformRandomNumberGenerator {
public:
    virtual unsigned operator()() = 0;
    virtual unsigned max() const = 0;
    virtual unsigned min() const = 0;
};

class SdtUNRG01 : public UniformRandomNumberGenerator {
    std::default_random_engine m_re;
public:
    SdtUNRG01(unsigned seed = 0): m_re(seed) {}

    unsigned operator()() override { return m_re(); }

    unsigned max() const override { return std::default_random_engine::max(); }
    unsigned min() const override { return std::default_random_engine::min(); }
};



class Noise {
    tvb::Float m_ntau;
    tvb::Float m_dt_sqrt;
    UniformRandomNumberGenerator* m_random_stream;
    std::normal_distribution<tvb::Float> m_normal_dist;

public:
    Noise(tvb::Float dt, tvb::Float ntau = 0.0, UniformRandomNumberGenerator* urng = NULL):
    m_ntau(ntau),
    m_normal_dist(0.0, 1.0)
    {
        m_dt_sqrt = sqrt(dt);
        if (urng == NULL)
            m_random_stream = new SdtUNRG01();
        else
            m_random_stream = urng;
    }

    tvb::TArray2d generate(int rows, int cols) {
        if (m_ntau > 0.0)
            // TODO: implement coloured noise
            throw("Coloured noise not implemented!");
        else
            return this->white(rows, cols);
    }

    tvb::TArray2d white(int rows, int cols) {
        tvb::TArray2d random_matrix = tvb::TArray2d::NullaryExpr(rows, cols, [this]() { return this->m_normal_dist(*this->m_random_stream); });
        return random_matrix * m_dt_sqrt;
    }

    virtual tvb::TArray1d gfun(const tvb::TArray2d& X) = 0;

};

class Additive : public Noise {
    tvb::TArray1d m_sqrt_2nsig;

public:
    Additive(const tvb::TArray1d& nsig, double dt, double ntau = 0.0): Noise(dt, ntau) {
        m_sqrt_2nsig = sqrt((2.0 * nsig).array());
    }

    tvb::TArray1d gfun(const tvb::TArray2d& X) override {
        return m_sqrt_2nsig;
    }

};


#endif //TVB_CPP_NOISE_H
