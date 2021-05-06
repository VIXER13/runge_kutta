#ifndef RUNGE_KUTTA_COEFFICIENTS_HPP
#define RUNGE_KUTTA_COEFFICIENTS_HPP

#include <array>

namespace ode {

template<class T>
struct triplet {
    size_t i, j;
    T val;
};

template<class T>
class runge_kutta4 {
    static_assert(std::is_floating_point<T>{}, "T must be floating point.");
    using p = std::pair<size_t, T>;
    using t = triplet<T>;

public:
    static constexpr uintmax_t order = 4;
    static constexpr size_t stages = 4;
    static constexpr bool nested = false;
    static constexpr std::array<t, 3> aij = {
        t{1, 0, T{1}/T{2}},
                            t{2, 1, T{1}/T{2}},
                                                t{3, 2,      T{1}}
    };
    static constexpr std::array<p, 4> bj = {
        p{0, T{1}/T{6}}, p{0, T{1}/T{6}}, p{0, T{1}/T{6}}, p{0, T{1}/T{6}},
    };
    static constexpr std::array<T, stages> ci = {T{0}, T{1}/T{2}, T{1}/T{2}, T{1}};
};

template<class T>
class dormand_prince {
    static_assert(std::is_floating_point<T>{}, "T must be floating point.");
    using p = std::pair<size_t, T>;
    using t = triplet<T>;

public:
    static constexpr uintmax_t order = 5;
    static constexpr size_t stages = 7;
    static constexpr bool nested = true;
    static constexpr std::array<t, 20> aij = {
        t{1, 0, T{    1}/T{   5}},
        t{2, 0, T{    3}/T{  40}}, t{2, 1, T{     9}/T{  40}},
        t{3, 0, T{   44}/T{  45}}, t{3, 1, T{   -56}/T{  15}}, t{3, 2, T{   32}/T{   9}},
        t{4, 0, T{19372}/T{6561}}, t{4, 1, T{-25360}/T{2187}}, t{4, 2, T{64448}/T{6561}}, t{4, 3, T{-212}/T{729}},
        t{5, 0, T{ 9017}/T{3168}}, t{5, 1, T{  -355}/T{  33}}, t{5, 2, T{46732}/T{5247}}, t{5, 3, T{  49}/T{176}}, t{5, 4, T{ -5103}/T{ 18656}},
        t{6, 0, T{   35}/T{ 384}},                             t{6, 2, T{  500}/T{1113}}, t{6, 3, T{ 125}/T{192}}, t{6, 4, T{ -2187}/T{  6784}}, t{6, 5, T{ 11}/T{  84}}
    };
    static constexpr std::array<p, 5> bj = {
        p{0, T{  35}/T{  384}}, p{2, T{ 500}/T{ 1113}}, p{3, T{ 125}/T{192}}, p{4, T{ -2187}/T{  6784}}, p{5, T{ 11}/T{  84}}
    };
    static constexpr std::array<p, 6> bj_ = {
        p{0, T{5179}/T{57600}}, p{2, T{7571}/T{16695}}, p{3, T{ 393}/T{640}}, p{4, T{-92097}/T{339200}}, p{5, T{187}/T{2100}}, p{6, T{1}/T{40}}
    };
    static constexpr std::array<T, stages> ci = {T{0}, T{1}/T{5}, T{3}/T{10}, T{4}/T{5}, T{8}/T{9}, T{1}, T{1}};
};

}

#endif