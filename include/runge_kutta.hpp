#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include "utils.hpp"
#include <limits>
#include <tuple>
#include <vector>

template<class T>
class runge_kutta4 {
    static_assert(std::is_floating_point<T>{}, "T must be floating point.");

    using p = std::pair<size_t, T>;
    template<size_t N>
    using a = std::array<p, N>;

public:
    static constexpr auto aij = std::make_tuple(
        a<1>{p{0,        T{0}}},
        a<1>{p{0, T{1} / T{2}}},
        a<1>{p{1, T{1} / T{2}}},
        a<1>{p{2,        T{1}}}
    );

    static constexpr size_t stages = std::tuple_size<decltype(aij)>{};

    static constexpr std::array<T, stages> ci = {T{0}, T{1} / T{2}, T{1} / T{2}, T{1}};
    static constexpr std::array<T, stages> bj = {T{1} / T{6}, T{1} / T{3}, T{1} / T{3}, T{1} / T{6}};
};

template<class T>
struct parameters {
    T tau_min  = 0,
      tau_init = T{1} / T{10},
      tau_max  = std::numeric_limits<T>::max(),
      tol      = T{1} / T{1000},
      fact_min = T{1} / T{2},
      fact     = T{1},
      fact_max = T{2};
    std::array<T, 2> time_interval = {T{0}, T{1}};
};

class _runge_kutta {
    template<size_t I, class T, template<class> class Coeffs, class Container, class System>
    static void stage(std::array<Container, Coeffs<T>::stages>& k,
                      const Container& sol, const T time, const T tau, const System& system) {
        std::fill(k[I].begin(), k[I].end(), T{0});
        using namespace utils;
        for(const auto [ind, val] : std::get<I>(Coeffs<T>::aij))
            k[I] += val * k[ind];
        k[I] = system(time + tau * Coeffs<T>::ci[I], sol + tau * k[I]);
    }

    template<class T, template<class> class Coeffs, class Container, class System, size_t ...I>
    static void make_step_impl(std::array<Container, Coeffs<T>::stages>& k,
                               const Container& sol, const T time, const T tau, const System& system, std::index_sequence<I...>) {
        (stage<I, T, Coeffs>(k, sol, time, tau, system), ...);
    }

    template<class T, template<class> class Coeffs, class Container, class System>
    static void make_step(std::array<Container, Coeffs<T>::stages>& k,
                          const Container& sol, const T time, const T tau, const System& system) {
        make_step_impl<T, Coeffs>(k, sol, time, tau, system, std::make_index_sequence<Coeffs<T>::stages>{});
    }

public:
    template<class T, template<class> class Coeffs, class Container, class System>
    friend std::pair<std::vector<T>, std::vector<Container>>
    runge_kutta(const parameters<T>& params, const Container& init, const System& system);
};

template<class T, template<class> class Coeffs, class Container, class System>
std::pair<std::vector<T>, std::vector<Container>>
runge_kutta(const parameters<T>& params, const Container& init, const System& system) {
    Container K = init;
    std::vector<T> time(1, params.time_interval.front());
    std::vector<Container> sol(1, init);
    std::array<Container, Coeffs<T>::stages> k;
    k.fill(init);
    const uintmax_t steps = (params.time_interval.back() - params.time_interval.front()) / params.tau_init;
    for(uintmax_t i = 0; i < steps; ++i) {
        _runge_kutta::template make_step<T, Coeffs>(k, sol.back(), time.back(), params.tau_init, system);
        std::fill(K.begin(), K.end(), T{0});
        using namespace utils;
        for(size_t i = 0; i < Coeffs<T>::stages; ++i)
            K += Coeffs<T>::bj[i] * k[i];
        time.push_back(time.back() + params.tau_init);
        sol.push_back(sol.back() + params.tau_init * K);
    }
    return {std::move(time), std::move(sol)};
}

#endif