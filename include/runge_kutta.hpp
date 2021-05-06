#ifndef RUNGE_KUTTA_RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_RUNGE_KUTTA_HPP

#include "coefficients.hpp"
#include "utils.hpp"
#include <vector>
#include <numeric>

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
    explicit _runge_kutta() = default;

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

    template<class T, template<class> class Coeffs, class Container>
    static void calc_sol_step(std::array<Container, 2>& K, const std::array<Container, Coeffs<T>::stages>& k) {
        using namespace utils;
        std::fill(K[0].begin(), K[0].end(), 0);
        for(const auto& [j, coeff] : Coeffs<T>::bj)
            K[0] += coeff * k[j];
        std::fill(K[1].begin(), K[1].end(), 0);
        for(const auto& [j, coeff] : Coeffs<T>::bj_)
            K[1] += coeff * k[j];
    }

    template<class T, template<class> class Coeffs, class Container>
    static void calc_sol_step(Container& K, const std::array<Container, Coeffs<T>::stages>& k) {
        using namespace utils;
        std::fill(K.begin(), K.end(), 0);
        for(const auto& [j, coeff] : Coeffs<T>::bj)
            K += coeff * k[j];
    }

    template<uintmax_t Order, bool Nested, class T, class Container>
    static T calc_error(std::array<Container, 2>& K) {
        T err = 0;
        for(size_t i = 0; i < K[0].size(); ++i)
            if(const T temp = std::abs((K[0][i] - K[1][i]) / std::max(K[0][i], K[1][i])); temp > err)
                err = temp;
        return Nested ? err : err / (utils::power<Order>(2) - 1);
    }

    template<class T, class Container>
    static void calc_new_sol(Container& K, const Container& sol, const T tau) {
        using namespace utils;
        K *= tau;
        K += sol;
    }

    template<uintmax_t Order, class T>
    static T tau_factor(const runge_kutta_parameters<T>& parameters, const T err) {
        return std::min(parameters.fact_max,
               std::max(parameters.fact_min,
                        parameters.fact * std::pow(parameters.tol / err, T{1} / T{Order + 1})));
    }

public:
    template<class T, template<class> class Coeffs, class Container, class System>
    friend std::tuple<std::vector<T>, std::vector<Container>, std::vector<T>, uintmax_t>
    runge_kutta(const runge_kutta_parameters<T>& parameters, const Container& init, const System& system);
};

template<class T, template<class> class Coeffs, class Container, class System>
std::tuple<std::vector<T>, std::vector<Container>, std::vector<T>, uintmax_t>
runge_kutta(const runge_kutta_parameters<T>& parameters, const Container& init, const System& system) {
    uintmax_t defect_step = 0;
    std::vector<T> loc_err(1, 0);
    T tau = parameters.tau_init;
    std::vector<Container> sol(1, init);
    std::vector<T> time(1, parameters.time_interval.front());
    std::array<Container, Coeffs<T>::stages> k; k.fill(init);
    std::array<Container, 2> K; K.fill(init);
    while (time.back() < parameters.time_interval.back()) {
        bool step_resize = parameters.autoresize;
        while (true) {
            _runge_kutta::template make_step<T, Coeffs>(k, time.back(), sol.back(), tau, system);
            if (!step_resize) {
                _runge_kutta::template calc_sol_step<T, Coeffs>(K[0], k);
                _runge_kutta::template calc_new_sol(K[0], sol.back(), tau);
                loc_err.push_back(0);
                time.push_back(time.back() + tau);
                sol.push_back(K[0]);
                break;
            }

            if constexpr (Coeffs<T>::nested) {
                _runge_kutta::template calc_sol_step<T, Coeffs>(K, k);
                _runge_kutta::template calc_new_sol(K[0], sol.back(), tau);
                _runge_kutta::template calc_new_sol(K[1], sol.back(), tau);
                const T err = _runge_kutta::template calc_error<Coeffs<T>::order, Coeffs<T>::nested, T>(K);
                if (err < parameters.tol) {
                    loc_err.push_back(err);
                    time.push_back(time.back() + tau);
                    sol.push_back(K[0]);
                }
                tau *= _runge_kutta::template tau_factor<Coeffs<T>::order>(parameters, err);
                if (err < parameters.tol)
                    break;
            } else {
                _runge_kutta::template calc_sol_step<T, Coeffs>(K[0], k);
                _runge_kutta::template calc_new_sol(K[0], sol.back(), tau);

                const T half_tau = 0.5 * tau;
                _runge_kutta::template make_step<T, Coeffs>(k, time.back(), sol.back(), half_tau, system);
                _runge_kutta::template calc_sol_step<T, Coeffs>(K[1], k);
                _runge_kutta::template calc_new_sol(K[1], sol.back(), half_tau);
                const Container half_sol = K[1];
                _runge_kutta::template make_step<T, Coeffs>(k, time.back() + half_tau, half_sol, half_tau, system);
                _runge_kutta::template calc_sol_step<T, Coeffs>(K[1], k);
                _runge_kutta::template calc_new_sol(K[1], half_sol, half_tau);

                const T err = _runge_kutta::template calc_error<Coeffs<T>::order, Coeffs<T>::nested, T>(K);
                if (err < parameters.tol_small)
                    tau *= parameters.fact_max;
                else if (err < parameters.tol) {
                    loc_err.push_back(err);
                    time.push_back(time.back() + tau);
                    sol.push_back(K[0]);
                    break;
                } else
                    tau *= parameters.fact_min;
            }

            ++defect_step;
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
    }
    return {std::move(time), std::move(sol), std::move(loc_err), defect_step};
}

}

#endif