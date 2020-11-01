#pragma once

#include <vector>

#include <iostream>

namespace simplex_method {

    template<typename T>
    using Matrix = std::vector<std::vector<T>>;

    template<typename T>
    Matrix<T> transposedCopy(Matrix<T> const& matrix) {
        auto const rows = matrix.size(); // > 0
        if (rows == 0) return {};

        auto const columns = matrix[0].size();
        Matrix<T> result(columns);
        for (size_t i = 0; i < columns; ++i) result[i] = std::vector<T>(rows);

        for (size_t i = 0; i < rows; ++i) {
            auto const& row = matrix[i];
            for (size_t j = 0; j < columns; ++j) result[j][i] = row[j];
        }

        return result;
    }
}
