#include <iostream>
#include <cmath>
#include "array_extended.hpp"
#include "dormand_prince.hpp"


int main() {
    auto [time, sol] = dormand_prince<double, 2>({0., 1.}, {0., 2.*M_PI}, 0.01, 
                        [](double t, std::array<double, 2>)->std::array<double, 2> { return {cos(t), -20.*sin(20.*t)}; },
                        {.tol = 1e-10});

    for(size_t i = 0; i < sol.size(); ++i)
        std::cout << "step = " << i << " time = " << time[i] << " sol = " << sol[i] << std::endl;
        
    return 0;
}