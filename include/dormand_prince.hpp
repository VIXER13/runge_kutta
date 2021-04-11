#ifndef DORMAND_PRINCE_HPP
#define DORMAND_PRINCE_HPP

#include "utils.hpp"
#include <vector>
#include <array>
#include <tuple>

template<size_t I, class T, class Container, class System, size_t Size>
void stage(Container& K, std::array<Container, 7>& k,
           const T time, const T tau, const Container& sol,
           const std::array<T, Size>& aij, const T c) {
    std::fill(K.begin(), K.end(), T{0});
    using namespace utils;
    for(size_t j = 0; j < I; ++j)
        K += aij[j] * k[j];
    k[I] = system(time + tau * c, sol + tau * K);
}

template<class T, class Container, class System, class ...Coeffs, size_t ...I>
void calc_k_impl(std::index_sequence<I...>) {
    (stage<I>(), ...);
}

template<class T, class Container, class System, class ...Coeffs>
void calc_k(const std::tuple<Coeffs...>& aij) {
    calc_k_impl(std::make_index_sequence<sizeof...(Coeffs)>());
}

template<class T, class Container, class System>
std::pair<std::vector<T>, std::vector<Container>>
dormand_prince(const Container& init, const std::array<T, 2>& time_interval, const T tau, const System& system) {
    static constexpr auto
        butcher_aij = std::make_tuple(
            std::array<T, 1>{T{0}},
            std::array<T, 1>{T{    1}/T{   5}},
            std::array<T, 2>{T{    3}/T{  40}, T{     9}/T{  40}},
            std::array<T, 3>{T{   44}/T{  45}, T{   -56}/T{  15}, T{   32}/T{   9}},
            std::array<T, 4>{T{19372}/T{6561}, T{-25360}/T{2187}, T{64448}/T{6561}, T{-212}/T{729}},
            std::array<T, 5>{T{ 9017}/T{3168}, T{  -355}/T{  33}, T{46732}/T{5247}, T{  49}/T{176}, T{-5103}/T{18656}},
            std::array<T, 6>{T{   35}/T{ 384},              T{0}, T{  500}/T{1113}, T{ 125}/T{192}, T{-2187}/T{ 6784}, T{11}/T{84}}
        );

    static constexpr std::array<T, 7>
        butcher_c  = {            T{0}, T{1}/T{5}, T{   3}/T{   10}, T{  4}/T{  5}, T{     8}/T{     9},           T{1},       T{1}},
        butcher_b  = {T{  35}/T{  384},      T{0}, T{ 500}/T{ 1113}, T{125}/T{192}, T{ -2187}/T{  6784}, T{ 11}/T{  84},       T{0}},
        butcher_b_ = {T{5179}/T{57600},      T{0}, T{7571}/T{16695}, T{393}/T{640}, T{-92097}/T{339200}, T{187}/T{2100}, T{1}/T{40}};

    auto k = std::array<Container, 7>{}.fill(init);

    std::vector<T> time(1, 0);
    std::vector<Container> sol(1, init);
    Container K_pred = init, K_corr = init;
    while(time.back() < time_interval[1]) {
        while(true) {
//            for(size_t i = 0; i < 7; ++i) {
//                std::fill(K_pred.begin(), K_pred.end(), 0);
//                for(size_t j = 0; j < i; ++j)
//                    K_pred += std::get<i>(butcher_table_coeff)[j] * k[j];
//                k[i] = system(time.back() + tau * butcher_table_tau[i], sol.back() + tau * K_pred);
//            }


            std::fill(K_pred.begin(), K_pred.end(), T{0});
            std::fill(K_corr.begin(), K_corr.end(), T{0});
            for(size_t i = 0; i < k.size(); ++i) {
                using namespace utils;
                K_pred += butcher_b [i] * k[i];
                K_corr += butcher_b_[i] * k[i];
            }
            break;
        }
    }

    return {std::move(time), std::move(sol)};
}

#endif