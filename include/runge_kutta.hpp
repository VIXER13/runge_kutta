#ifndef RUNGE_KUTTA_RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_RUNGE_KUTTA_HPP

#include "coefficients.hpp"
#include "utils.hpp"
#include <vector>
#include <numeric>
#include <iostream>

namespace ode {

template<class T>
struct runge_kutta_parameters {
    T tau_min   = 0,                               // Минимальный размер шага
      tau_init  = T{1} / T{10},                    // Начальный размер шага
      tau_max   = std::numeric_limits<T>::max(),   // Максимальный размер шага
      tol       = T{1} / T{1000},                  // Точность
      tol_small = T{1} / T{1000000},               //
      fact_min  = T{1} / T{2},                     // Минимальный фактор
      fact      = T{1},                            // Гарантийный фактор
      fact_max  = T{2};                            // Максимальный фактор
    std::array<T, 2> time_interval = {T{0}, T{1}}; // Интервал интегрирования
};

class _runge_kutta {
    template<class T, template<class> class Coeffs, class Container, class System>
    static void make_step(std::array<Container, Coeffs<T>::stages>& k,
                          const T time, const Container& sol, const T tau, const System& system) {
        size_t curr = 0;
        using namespace utils;
        for(const auto& [i, j, val] : Coeffs<T>::aij) {
            if (i == curr)
                k[i] += val * k[j];
            else {
                k[curr] = system(time + tau * Coeffs<T>::ci[curr], sol + tau * k[curr]);
                k[i] = val * k[j];
                curr = i;
            }
        }
        k[curr] = system(time + tau * Coeffs<T>::ci[curr], sol + tau * k[curr]);
    }

    template<class T, template<class> class Coeffs, size_t Index, class Container>
    static void calc_sol_step(std::array<Container, 2>& K,
                              const std::array<Container, Coeffs<T>::stages>& k) {
        size_t curr = -1;
        for(const auto& [i, j, val] : Coeffs<T>::bj) {
            using namespace utils;
            if (i == curr)
                K[Coeffs<T>::nested ? i : Index] += val * k[j];
            else {
                K[Coeffs<T>::nested ? i : Index] = val * k[j];
                curr = i;
            }
        }
    }

    template<uintmax_t Order, bool Nested, class T, class Container>
    static T calc_error(std::array<Container, 2>& K) {
        T err = 0;
        for(size_t i = 0; i < K[0].size(); ++i) {
            const T temp = Nested ? std::abs((K[0][i] - K[1][i]) / K[1][i]) : std::abs(K[0][i] - K[1][i]);
            if(temp > err)
                err = temp;
        }
        static constexpr T factor = T{1} / (utils::power<Order>(2) - 1);
        return factor * err;
    }

    template<uintmax_t Order, class T>
    static T tau_factor(const runge_kutta_parameters<T>& parameters, const T err) {
        static constexpr T power = T{1} / T{Order + 1};
        return std::min(parameters.fact_max,
               std::max(parameters.fact_min,
                        parameters.fact * std::pow(parameters.tol / err, power)));
    }

public:
    template<class T, template<class> class Coeffs, class Container, class System>
    friend std::pair<std::vector<T>, std::vector<Container>>
    runge_kutta(const runge_kutta_parameters<T>& parameters, const Container& init, const System& system);
};

template<class T, template<class> class Coeffs, class Container, class System>
std::pair<std::vector<T>, std::vector<Container>>
runge_kutta(const runge_kutta_parameters<T>& parameters, const Container& init, const System& system) {
    T tau = parameters.tau_init;
    std::vector<Container> sol(1, init);
    std::vector<T> time(1, parameters.time_interval.front());
    std::array<Container, Coeffs<T>::stages> k; k.fill(init);
    std::array<Container, 2> K; K.fill(init);
    while (time.back() < parameters.time_interval.back()) {
        bool step_resize = true;
        while (true) {
            _runge_kutta::template make_step<T, Coeffs>(k, time.back(), sol.back(), tau, system);
            _runge_kutta::template calc_sol_step<T, Coeffs, 0>(K, k);

            if (!step_resize)
                break;

            if constexpr (Coeffs<T>::nested) {
                const T err = _runge_kutta::template calc_error<Coeffs<T>::order, Coeffs<T>::nested, T>(K);
                if (err < parameters.tol)
                    break;
                tau *= _runge_kutta::template tau_factor<Coeffs<T>::order>(parameters, err);
            } else {
                _runge_kutta::template make_step<T, Coeffs>(k, time.back(), sol.back(), 0.5 * tau, system);
                _runge_kutta::template calc_sol_step<T, Coeffs, 1>(K, k);
                const T err = _runge_kutta::template calc_error<Coeffs<T>::order, Coeffs<T>::nested, T>(K);
                if (err < parameters.tol_min)
                    tau *= parameters.fact_max;
                else if (err < parameters.tol)
                    break;
                tau *= parameters.fact_min;
            }

            if (tau < parameters.tau_min) {
                tau = parameters.tau_min;
                step_resize = false;
            }
            if (tau > parameters.tau_max) {
                tau = parameters.tau_max;
                step_resize = false;
            }
            if (time.back() + tau > parameters.time_interval.back()) {
                tau = parameters.time_interval.back() - time.back();
                step_resize = false;
            }
        }
        using namespace utils;
        time.push_back(time.back() + tau);
        sol.emplace_back(sol.back() + tau * K[0]);
    }
    return {std::move(time), std::move(sol)};
}

}

#endif