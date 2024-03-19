//
// Created by natxm on 19/03/2024.
//

#include "factory.h"

#include <tvb-cpp/simulator/models/reduced_ww_ext.h>
#include <tvb-cpp/simulator/models/zerlaut.h>
#include <tvb-cpp/simulator/models/montbrio.h>

tvb::Model *tvb::Factory::new_model(const std::string &model_name, int N) {
    if (model_name == "ReducedWongWangExcInh") {
        return new ReducedWongWangExcInh(N);
    }
    else if (model_name == "ZerlautAdaptationFirstOrder") {
        return new ZerlautAdaptationFirstOrder(N);
    }
    else if (model_name == "ZerlautAdaptationSecondOrder") {
        return new ZerlautAdaptationSecondOrder(N);
    }
    else if (model_name == "Montbrio") {
        return new Montbrio(N);
    }
    throw std::runtime_error(string_format("Factory error: model name <%s> not available", model_name.c_str()));
}
