#pragma once

#include <iostream>
#include <simplex_method.hpp>
#include <simplex_method_matrix.hpp>
#include <stdexcept>
#include <vector>

namespace simplex_method {

    // F = cx -> max
    // Ax <= b
    // x >= 0
    template<typename T>
    [[nodiscard]] NormalizedOptimizationTask<T> normalizedOptimizationTask(bool max, std::vector<T> const& c,
                                                                           std::vector<std::vector<T>> const& a,
                                                                           std::vector<T> const& b) {
        auto const freeVariables = c.size();// n - m
        if (freeVariables == 0) throw std::invalid_argument("c cannot be empty");
        if (a.size() != freeVariables) throw std::invalid_argument("a and c should all be of the same length");

        auto const m = a[0].size();

        Matrix<T> conditions(m);
        {
            auto const rowLength = freeVariables + 2;
            size_t i = 0;
            for (auto const& line : a) {
                auto &condition = conditions[i];

                condition.reserve(rowLength);
                for (auto const& coefficient : line) condition.push_back(coefficient);
                condition.push_back(1);// basis variable
                condition.push_back(b[i++]);
            }
        }

        return NormalizedOptimizationTask(max, freeVariables + m, m, c, conditions);
    }
}// namespace simplex_method
