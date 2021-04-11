#ifndef RUNGE_KUTTA_RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_RUNGE_KUTTA_HPP

#include "coefficients.hpp"
#include "utils.hpp"
#include <vector>
#include <iostream>

namespace ode {

template<class T>
struct runge_kutta_parameters {
    T tau_min  = 0,                                // Минимальный размер шага
    tau_init = T{1} / T{10},                     // Начальный размер шага
    tau_max  = std::numeric_limits<T>::max(),    // Максимальный размер шага
    tol      = T{1} / T{1000},                   // Точность
    fact_min = T{1} / T{2},                      // Минимальный фактор
    fact     = T{1},                             // Гарантийный фактор
    fact_max = T{2};                             // Максимальный фактор
    std::array<T, 2> time_interval = {T{0}, T{1}}; // Интервал интегрирования
};

template<class T, template<class> class Coeffs, class Container, class System>
std::pair<std::vector<T>, std::vector<Container>>
runge_kutta(const runge_kutta_parameters<T>& parameters, const Container& init, const System& system) {
    std::vector<Container> sol(1, init);
    std::vector<T> time(1, parameters.time_interval.front());
    std::array<Container, Coeffs<T>::stages> k; k.fill(init);
    std::array<Container, Coeffs<T>::nested ? 2 : 1> K; K.fill(init);
    const T tau = parameters.tau_init;
    while (time.back() < parameters.time_interval.back()) {
        using namespace utils;
        size_t curr = 0;
        for(const auto& [i, j, val] : Coeffs<T>::aij) {
            if (i == curr) {
                k[i] += val * k[j];
            } else {
                k[curr] = system(time.back() + tau * Coeffs<T>::ci[curr], sol.back() + tau * k[curr]);
                k[i] = val * k[j];
                curr = i;
            }
        }
        k.back() = system(time.back() + tau * Coeffs<T>::ci[curr], sol.back() + tau * k.back());

        curr = -1;
        for(const auto& [i, j, val] : Coeffs<T>::bj) {
            if (i == curr)
                K[i] += val * k[j];
            else {
                K[i] = val * k[j];
                curr = i;
            }
        }

        time.push_back(time.back() + parameters.tau_init);
        sol.push_back(sol.back() + parameters.tau_init * K.front());
    }
    return {std::move(time), std::move(sol)};
}

}

#endif