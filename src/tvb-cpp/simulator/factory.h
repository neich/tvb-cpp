//
// Created by natxm on 19/03/2024.
//

#ifndef TVB_CPP_EXAMPLES_FACTORY_H
#define TVB_CPP_EXAMPLES_FACTORY_H

#include "model.h"

namespace tvb {

    class Factory {
    public:
        static Model* new_model(const std::string& model_name, int N);
    };
}


#endif //TVB_CPP_EXAMPLES_FACTORY_H
