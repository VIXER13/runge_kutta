#ifndef RUNGE_KUTTA4_HPP
#define RUNGE_KUTTA4_HPP

#include <utility>
#include "utils.hpp"

template<class T, class Container, class System>
std::pair<std::vector<T>, std::vector<Container>>
runge_kutta4(const Container& init, const std::array<T, 2>& time_interval, const T tau, const System& system) {
    Container K, k;
    const T tau2 = tau / 2, tau6 = tau / 6;
    const uintmax_t steps = (time_interval[1] - time_interval[0]) / tau;
    std::vector<T>         time(steps + 1); time[0] = time_interval[0];
    std::vector<Container> sol (steps + 1); sol [0]  = init;
    for(uintmax_t i = 0; i < steps; ++i) {
        using namespace utils;
        k = system(time[i],        sol[i]           ); K  = k;     // k1
        k = system(time[i] + tau2, sol[i] + tau2 * k); K += 2 * k; // k2
        k = system(time[i] + tau2, sol[i] + tau2 * k); K += 2 * k; // k3
        k = system(time[i] + tau,  sol[i] + tau  * k); K += k;     // k4
        time[i + 1] = time[i] + tau;
        sol [i + 1] = sol [i] + tau6 * K;
    }
    return {std::move(time), std::move(sol)};
}

#endif