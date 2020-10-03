#pragma once

#include <iostream>
#include <simplex_method_definitions.hpp>
#include <variant>
#include <vector>

namespace simplex_method {

    template<typename T>
    struct SuccessfulOptimizationResult {
        T value;
        std::vector<T> x;
    };

    template<typename T>
    using OptimizationResult = std::variant<SuccessfulOptimizationResult<T>, std::string>;
}// namespace simplex_method

template<typename T>
std::ostream &operator<<(std::ostream &out, simplex_method::SuccessfulOptimizationResult<T> result) {
    out << "F = " << result.value << "; {";

    auto const size = result.x.size();
    {
        size_t i = 0;
        for (auto const &value : result.x) {
            out << "x" << ++i << " = " << value;
            if (i != size) out << ", ";
        }
    }

    return out << "}";
}

template<typename T>
std::ostream &operator<<(std::ostream &out, simplex_method::OptimizationResult<T> result) {
    return std::holds_alternative<std::string>(result)
                           ? out << "[no result]: " << std::get<std::string>(result)
                           : out << std::get<simplex_method::SuccessfulOptimizationResult<T>>(result);
}
