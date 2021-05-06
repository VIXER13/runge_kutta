#ifndef RUNGE_KUTTA_UTILS_HPP
#define RUNGE_KUTTA_UTILS_HPP

#include <cstddef>
#include <tuple>
#include <cmath>
#include <utility>

namespace utils {

template<intmax_t N, class T>
constexpr T power(const T x) {
    if constexpr (N > 0 && N % 2 == 0) {
        const T temp = power<N / 2>(x);
        return temp * temp;
    } else if constexpr (N > 0 && N % 2 == 1)
        return x * power<N - 1>(x);
    else if constexpr (N < 0)
        return 1 / power<-N>(x);
    return 1;
}

template<class T, class Container>
Container operator*(const T& val, const Container& data) {
    Container new_data = data;
    for(size_t i = 0; i < data.size(); ++i)
        new_data[i] *= val;
    return std::move(new_data);
}

template<class Container>
Container operator+(const Container& lhs, const Container& rhs) {
    Container new_arr = lhs;
    for(size_t i = 0; i < new_arr.size(); ++i)
        new_arr[i] += rhs[i];
    return std::move(new_arr);
}

template<class Container>
Container& operator+=(Container& lhs, const Container& rhs) {
    for(size_t i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];
    return lhs;
}

template<class Container, class T>
Container& operator*=(Container& lhs, const T& val) {
    for(auto& it : lhs)
        it *= val;
    return lhs;
}

}

#endif