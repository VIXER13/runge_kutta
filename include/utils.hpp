#ifndef RUNGE_KUTTA_UTILS_HPP
#define RUNGE_KUTTA_UTILS_HPP

namespace utils {

template<class T, class Container>
Container operator*(const T val, const Container& arr) {
    Container new_arr = arr;
    for(size_t i = 0; i < arr.size(); ++i)
        new_arr[i] *= val;
    return std::move(new_arr);
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

}

#endif