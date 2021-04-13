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
    bool autoresize = true;
};

class _runge_kutta {
    enum class STEP_RESIZE_FLAG : uint8_t {FIX, RESIZE, NO};

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
            const T temp = std::abs(Nested ? (K[0][i] - K[1][i]) / std::max(K[0][i], K[1][i]) :
                                              K[0][i] - K[1][i]);
            if(temp > err)
                err = temp;
        }
        static constexpr T factor = T{1} / (utils::power<Order>(2) - 1);
        return factor * err;
    }

    template<class T, class Container>
    static void calc_new_sol(std::array<Container, 2>& new_sol, const Container& sol,
                             const std::array<Container, 2>& K, const T tau) {
        using namespace utils;
        for(size_t i = 0; i < 2; ++i)
            new_sol[i] = sol + tau * K[i];
    }

    template<uintmax_t Order, class T>
    static T tau_factor(const runge_kutta_parameters<T>& parameters, const T err) {
        static constexpr T power = T{1} / T{Order + 1};
        return std::min(parameters.fact_max,
               std::max(parameters.fact_min,
                        parameters.fact * std::pow(parameters.tol / err, power)));
    }

    template<class T, template<class> class Coeffs, class Container, class System>
    static STEP_RESIZE_FLAG resize_tau(T& tau, std::array<Container, 2>& new_sol,
                                       const runge_kutta_parameters<T>& parameters,
                                       const T time, const Container& x,
                                       const std::array<Container, 2>& K,
                                       const std::array<Container, Coeffs<T>::stages>& k) {

    }

public:
    template<class T, template<class> class Coeffs, class Container, class System>
    friend std::tuple<std::vector<T>, std::vector<Container>, std::vector<T>>
    runge_kutta(const runge_kutta_parameters<T>& parameters, const Container& init, const System& system);
};

template<class T, template<class> class Coeffs, class Container, class System>
std::tuple<std::vector<T>, std::vector<Container>, std::vector<T>>
runge_kutta(const runge_kutta_parameters<T>& parameters, const Container& init, const System& system) {
    T tau = parameters.tau_init;
    std::vector<Container> sol(1, init);
    std::vector<T> time(1, parameters.time_interval.front());
    std::array<Container, Coeffs<T>::stages> k; k.fill(init);
    std::array<Container, 2> K; K.fill(init);
    std::array<Container, 2> new_sol;
    std::vector<T> loc_err(1, 0);
    while (time.back() < parameters.time_interval.back()) {
        bool step_resize = parameters.autoresize;
        while (true) {
            _runge_kutta::template make_step<T, Coeffs>(k, time.back(), sol.back(), tau, system);
            _runge_kutta::template calc_sol_step<T, Coeffs, 0>(K, k);

            if (!step_resize) {
                _runge_kutta::template calc_new_sol(new_sol, sol.back(), K, tau);
                break;
            }

            if constexpr (Coeffs<T>::nested) {
                _runge_kutta::template calc_new_sol(new_sol, sol.back(), K, tau);
                const T err = _runge_kutta::template calc_error<Coeffs<T>::order, Coeffs<T>::nested, T>(new_sol);
                if (err < parameters.tol)
                    break;
                tau *= _runge_kutta::template tau_factor<Coeffs<T>::order>(parameters, err);
            } else {
                _runge_kutta::template make_step<T, Coeffs>(k, time.back(), sol.back(), 0.5 * tau, system);
                _runge_kutta::template calc_sol_step<T, Coeffs, 1>(K, k);
                _runge_kutta::template calc_new_sol(new_sol, sol.back(), K, tau);
                const T err = _runge_kutta::template calc_error<Coeffs<T>::order, Coeffs<T>::nested, T>(new_sol);
                if (err < parameters.tol_small)
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
        loc_err.push_back(_runge_kutta::template calc_error<Coeffs<T>::order, Coeffs<T>::nested, T>(new_sol));
        time.push_back(time.back() + tau);
        sol.emplace_back(std::move(new_sol[0]));
    }
    return {std::move(time), std::move(sol), std::move(loc_err)};
}

}

#endif