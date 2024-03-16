//
// Created by natxm on 13/12/2023.
//

#ifndef TVB_CPP_PARAM_OBJECT_H
#define TVB_CPP_PARAM_OBJECT_H

#include <string>
#include <vector>

#include <tvb-cpp/definitions.h>


namespace tvb {

    class ParamObject {
    public:
        [[nodiscard]] virtual std::vector<std::string> get_param_list() const = 0;

        virtual void set_param(const std::string &param, tvb::Float value) = 0;

        [[nodiscard]] virtual tvb::Float get_param(const std::string &param) const = 0;

        virtual void init_dependant() {}

    };

}
#endif //TVB_CPP_PARAM_OBJECT_H
