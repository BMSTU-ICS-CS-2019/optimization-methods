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
// _________________________________________
// | 1 3 5 6 1 3 2 1 1 1 2 0 0 0.5 1 5 3 8 |
// |---------------------------------------|
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

    std::cout << "Enter `A`-vector: \n";
    simplex_method::Matrix<N> a(m);
    for (size_t i = 0; i < m; ++i) {
        (std::cout << "\tEnter line #" << i + 1 << ": ").flush();
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


    std::cout << "\n========= Direct task =========" << std::endl;

    auto const directTask = simplex_method::normalizedOptimizationTask(true, c, a, b);
    std::cout << "\t< Task >\n" << directTask << "\n" << std::endl;
    auto const directResult = directTask.compute(std::cout);
    std::cout << "\n\t< Result >\n" << directResult << "\n" << std::endl;

    std::cout << "===============================" << std::endl;


    std::cout << "\n========== Dual task ==========" << std::endl;

    auto const dualTask = simplex_method::normalizedDualOptimizationTask(true, c, a, b);
    std::cout << "\t< Task >\n" << dualTask << "\n" << std::endl;
    auto const dualResult = dualTask.compute(std::cout);
    std::cout << "\n\t< Result >\n" << dualResult << "\n" << std::endl;

    std::cout << "===============================" << std::endl;

    return !std::holds_alternative<simplex_method::SuccessfulOptimizationResult<N>>(directResult);
}
