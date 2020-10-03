#include <iostream>

#include <boost/rational.hpp>
#include <optional>
#include <simplex_method.hpp>
#include <simplex_task_normalizer.hpp>

template<typename T>
int performTask();

// F = cx -> max
// Ax <= b
// x >= 0
int main() {
    (std::cout << "Enter number-mode (1 for decimals, 2 for rationals): ").flush();

    char mode;
    std::cin >> mode;
    switch (mode) {
        case '1': return performTask<double>();
        case '2': return performTask<boost::rational<std::intmax_t>>();
        default: {
            std::cout << "Unsupported mode: " << mode << std::endl;

            return 1;
        }
    }
}

template<typename N>
int performTask() {
    (std::cout << "Enter the number of free variables (n - m): ").flush();
    size_t freeVariables;
    std::cin >> freeVariables;

    (std::cout << "Enter `c`-vector: ").flush();
    std::vector<N> c;
    c.reserve(freeVariables);
    for (size_t i = 0; i < freeVariables; ++i) {
        N value;
        std::cin >> value;
        c.push_back(value);
    }

    (std::cout << "Enter the number of conditions: ").flush();
    size_t m;
    std::cin >> m;

    (std::cout << "Enter `a`-vector: ").flush();
    simplex_method::Matrix<N> a(m);
    for (size_t i = 0; i < m; ++i) {
        auto &condition = a[i];
        for (size_t j = 0; j < freeVariables; ++j) {
            N value;
            std::cin >> value;
            condition.push_back(value);
        }
    }

    (std::cout << "Enter `b`-vector: ").flush();
    std::vector<N> b;
    b.reserve(freeVariables);
    for (size_t i = 0; i < m; ++i) {
        N value;
        std::cin >> value;
        b.push_back(value);
    }

    auto const task = simplex_method::normalizedOptimizationTask(true, c, a, b);
    auto const result = task.compute();

    const bool failure = !std::holds_alternative<simplex_method::SuccessfulOptimizationResult<N>>(result);
    std::cout << "Result:\n" << result << std::endl;

    return failure;
    /*
    auto normalized = simplex_method::normalizedOptimizationTask<double>(
            true, {2, 5, 3}, {{2, 1, 2}, {1, 2, 0}, {0, 0.5, 1}}, {6, 6, 2});
    std::cout << normalized.compute() << std::endl;
    std::cout << "Hello, World!"
              << simplex_method::NormalizedOptimizationTask<double>(false, 5, 3, {1, -1},
                                                                    {
                                                                            {1, -2, 1, 2},
                                                                            {2, -1, -1, 2},
                                                                            {1, 1, 1, 5},
                                                                    })
                         .compute()
              << std::endl;
              */

    return 0;
}


