#pragma once

#include <cmath>
#include <iostream>
#include <memory>
#include <optimization_result.hpp>
#include <optional>
#include <type_traits>
#include <variant>
#include <vector>

namespace simplex_method {

    template<typename T>
    class NormalizedOptimizationTask final {
        static_assert(std::is_default_constructible<T>::value, "T should be default-constructible");

        T ZERO = T();

        bool max_;
        std::size_t n_, m_;
        // F = coefficients[0] + sum({i = 1..n} coefficients[i] * x_i) -> max
        std::vector<T> coefficients_;
        // {i = 0..(m-1)} sum({j = 0..(n-m-1)} conditions[i][j] * x_j) + conditions[i][n - m] = conditions[n - m + 1]
        Matrix<T> conditions_;

    public:
        NormalizedOptimizationTask(bool max, size_t n, size_t m, const std::vector<T> &coefficients,
                                   const std::vector<std::vector<T>> &conditions);

        [[nodiscard]] OptimizationResult<T> compute(std::ostream &debugOut) const;

        std::ostream& appendTo(std::ostream& out) const;

    private:// helper methods
        void swapMatrixLines_(Matrix<T> &simplexMatrix, size_t iX, size_t jX) const;

        void appendVariableIndices_(std::vector<size_t> const& indices, std::ostream &out) const;
    };

    template<typename T>
    NormalizedOptimizationTask<T>::NormalizedOptimizationTask(bool max, size_t n, size_t m,
                                                              std::vector<T> const& coefficients,
                                                              std::vector<std::vector<T>> const& conditions)
        : max_(max), n_(n), m_(m), coefficients_(coefficients), conditions_(conditions) {}

    template<typename T>
    [[nodiscard]] OptimizationResult<T> NormalizedOptimizationTask<T>::compute(std::ostream &debugOut) const {
        // step 1: fill the initial matrix
        std::vector<size_t> indices(n_);
        for (size_t index = 0; index < n_; index++) indices[index] = index;

        std::vector<std::vector<T>> matrix(m_ + 1u);
        auto const freeVariables = n_ - m_, rowLength = freeVariables + 1u;
        for (size_t i = 0; i < m_; ++i) {
            auto& row = matrix[i];
            row.reserve(rowLength);

            auto const& conditions = conditions_[i];
            auto const divisor = conditions[freeVariables];
            row.push_back(conditions[rowLength] / divisor);
            for (size_t j = 0; j < freeVariables; ++j) row.push_back(conditions[j] / divisor);
        }
        {
            auto& row = matrix[m_];
            row.reserve(rowLength);
            row.push_back(ZERO); // zero
            for (auto const& coefficient : coefficients_) row.push_back(coefficient);
        }
        appendVariableIndices_(indices, debugOut << "Generated initial simplex matrix:\n" << matrix << "\n");
        debugOut << std::endl;

        // step 2: reach target solution
        do {
            std::optional<size_t> negativeRowIndex;
            for (size_t rowIndex = 0u; rowIndex < m_; ++rowIndex) {
                auto& row = matrix[rowIndex];
                auto const si0 = matrix[rowIndex][0];

                if (si0 < ZERO) {
                    negativeRowIndex = rowIndex;
                    break;
                }
            }
            if (!negativeRowIndex) break; // no more negative si0s

            auto &row = matrix[*negativeRowIndex];
            // find solving coordinates
            size_t resolvingColumn = 0, resolvingRow;
            for (size_t j = 1u; j < rowLength; ++j) if (row[j] < ZERO) {
                resolvingColumn = j;
                resolvingRow = 0;

                std::optional<T> ratio;
                for (size_t i = 0u; i < m_; ++i) {
                    auto const currentRatio = matrix[i][0] / matrix[i][resolvingColumn];
                    if (currentRatio > ZERO && (!ratio || currentRatio < *ratio)) {
                        ratio = currentRatio;
                        resolvingRow = i;
                    }
                }
                if (!ratio) return "Resolving row cannot be found";

                break; // resolving coordinates got found
            }
            if (resolvingColumn == 0) return "Resolving column cannot be found";

            debugOut << "Swapping row " << resolvingRow << " with column " << resolvingColumn << std::endl;
            std::swap(indices[resolvingColumn - 1], indices[rowLength - 1 + resolvingRow]);
            swapMatrixLines_(matrix, resolvingRow, resolvingColumn);

            appendVariableIndices_(indices, debugOut << "Resulting matrix:\n" << matrix << "\n");
            debugOut << std::endl;
        } while (true);
        debugOut << "Initial target solution has been found" << std::endl;

        // step 3: find the most appropriate value of F
        {
            auto const& lastRow = matrix[m_];

            bool doContinue;
            do {
                doContinue = false;
                for (size_t resolvingColumn = 1u; resolvingColumn < rowLength; ++resolvingColumn) if (
                    lastRow[resolvingColumn] > ZERO) {
                    size_t resolvingRow;

                    std::optional<T> ratio;
                    for (size_t i = 0u; i < m_; ++i) {
                        auto const currentRatio = matrix[i][0] / matrix[i][resolvingColumn];
                        if (currentRatio > ZERO && (!ratio || currentRatio < *ratio)) {
                            ratio = currentRatio;
                            resolvingRow = i;
                        }
                    }
                    if (!ratio) return "Resolving row cannot be found";

                    debugOut << "Swapping row " << resolvingRow << " with column " << resolvingColumn << std::endl;
                    std::swap(indices[resolvingColumn - 1], indices[rowLength - 1 + resolvingRow]);
                    swapMatrixLines_(matrix, resolvingRow, resolvingColumn);

                    appendVariableIndices_(indices, debugOut << "Resulting matrix:\n" << matrix << "\n");
                    debugOut << std::endl;

                    // current super-iteration has reached its end
                    doContinue = true;
                    break;
                }
            } while (doContinue);
        }
        debugOut << "Target solution has been found" << std::endl;

        std::vector<T> x(n_); // all values will be default-constructed, i.e. zeroed
        // for (size_t slot = 0; slot < freeVariables; ++slot) x[indices[slot]] = ZERO;
        for (size_t slot = freeVariables, i = 0; slot < n_; ++slot) { x[indices[slot]] = matrix.at(i++).at(0); }

        return SuccessfulOptimizationResult<T>{max_ ? -matrix[m_][0] : matrix[m_][0], x};
    }

    template <typename T>
    void NormalizedOptimizationTask<T>::swapMatrixLines_(Matrix<T>& simplexMatrix, const size_t iX,
                                                         const size_t jX) const {
        auto const rowLength = n_ - m_ + 2u;

        // center of rotation
        auto const& center = simplexMatrix[iX][jX];
        auto const invertedCenter = -center;

        // rotated row
        auto& rowIX = simplexMatrix[iX];
        //@formatter:off
        for (size_t i = 0; i <= m_; ++i) if (i != iX) {
            auto &row = simplexMatrix[i];
            for (size_t j = 0; j < rowLength; ++j) if (j != jX) row[j] -= rowIX[j] * row[jX] / center;
        }
        for (size_t j = 0; j < rowLength; ++j) if (j != jX) rowIX[j] /= center;
        for (size_t i = 0; i <= m_; ++i) if (i != iX) simplexMatrix[i][jX] /= invertedCenter;
        //@formatter:on

        rowIX[jX] = 1 / rowIX[jX];
    }

    template<typename T>
    void NormalizedOptimizationTask<T>::appendVariableIndices_(std::vector<size_t> const& indices,
                                                               std::ostream &out) const {
        out << "Free variables: {";
        size_t i = 0;
        {
            auto const freeVariables = n_ - m_;
            while (true) {
                out << "x" << indices[i] + 1;
                if (++i == freeVariables) break;
                else out << ", ";
            }
        }

        out << "}, Basis variables: {";
        while (true) {
            out << "x" << indices[i] + 1;
            if (++i == n_) break;
            else out << ", ";
        }
        out << "}";
    }

    template<typename T, typename P, typename S>
    inline static std::ostream& appendVariable_(std::ostream& out, T const& value, size_t const index,
                                                P prefix = "", S suffix = "") {
        if (value == 0) return out;

        out << prefix;
        if (value != 1) out << value;

        return out << 'x' << index << suffix;
    }

    template<typename T>
    inline static std::ostream& appendVariable_(std::ostream& out, T const& value, size_t const index) {
        return appendVariable_(out, value, index, "", "");
    }

    template <typename T>
    std::ostream& NormalizedOptimizationTask<T>::appendTo(std::ostream& out) const {
        // F = coefficients[0] + sum({i = 1..n} coefficients[i] * x_i) -> max
        // {i = 0..(m-1)} sum({j = 0..(n-m-1)} conditions[i][j] * x_j) + conditions[i][n - m] = conditions[n - m + 1]
        out << "F = ";

        {
            auto const size = coefficients_.size();
            for (size_t i = 0; i < size; ++i) {
                if (i != 0) out << " + ";
                appendVariable_(out, coefficients_[i], i + 1);
            }
        }
        out << " -> " << (this->max_ ? "max" : "min");

        {
            size_t basisIndex = n_ - m_;
            for (auto const& condition : this->conditions_) {
                out << "\n{ ";

                {
                    auto const preLastIndex = condition.size() - 2;
                    for (size_t i = 0; i < preLastIndex; ++i) appendVariable_(out, condition[i], i + 1, "", " + ");
                    appendVariable_(out, condition[preLastIndex], ++basisIndex) << " = " << condition[preLastIndex + 1];
                }
            }
        }

        return out;
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &out, Matrix<T> matrix) {
        auto const height = matrix.size();
        if (height == 0) return out << "{}";

        out << "{";
        for (auto const& row : matrix) {
            out << "\n\t";
            {
                size_t width = row.size();
                size_t j = 0;
                for (auto const& element : row) {
                    out << element;
                    if (++j != width) out << " ; ";
                }
            }
            //out << ((++i != height) ? "\n\t" : "\n");
        }

        return out << "\n}";
    }
}// namespace simplex_method

template<typename T>
std::ostream& operator<<(std::ostream& out, simplex_method::NormalizedOptimizationTask<T> const& task) {
    return task.appendTo(out);
}
