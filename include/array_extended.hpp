// В данном модуле добавлен дополнительный функционал для std::array.
// Данный функциона часто бывает удобен для отладки или некоторых мелких нужд.
// Он не претендует на полноту и оптимальность.

#ifndef ARRAY_EXTENDED_H
#define ARRAY_EXTENDED_H

#include <iostream>
#include <array>

template<class Type, size_t Size>
std::ostream& operator<<(std::ostream &os, const std::array<Type, Size> &vec)
{
    for(size_t i = 0; i < Size; ++i)
        os << vec[i] << " ";
    return os;
}

template<class Type, size_t Size>
std::istream& operator>>(std::istream &is, std::array<Type, Size> &vec)
{
    for(size_t i = 0; i < Size; ++i)
        is >> vec[i];
    return is;
}

template<class LeftType, class RightType, size_t Size>
std::array<decltype(LeftType() + RightType()), Size>
    operator+(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    std::array<decltype(LeftType() + RightType()), Size> c;
    for (size_t i = 0; i < Size; ++i)
        c[i] = a[i] + b[i];
    return c;
}

template<class LeftType, class RightType, size_t Size>
std::array<decltype(LeftType() - RightType()), Size>
    operator-(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    std::array<decltype(LeftType() - RightType()), Size> c;
    for (size_t i = 0; i < Size; ++i)
        c[i] = a[i] - b[i];
    return c;
}

template<class LeftType, class RightType, size_t Size>
decltype(LeftType() * RightType())
    operator*(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    decltype(LeftType() * RightType()) sum = decltype(LeftType() * RightType())();
    for(size_t i = 0; i < Size; ++i)
        sum += a[i] * b[i];
    return sum;
}

template<class ValType, class VecType, size_t Size>
std::array<decltype(VecType() * ValType()), Size>
    operator*(const std::array<VecType, Size> &a, const ValType &val)
{
    std::array<decltype(VecType() * ValType()), Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = a[i] * val;
    return c;
}

template<class ValType, class VecType, size_t Size>
std::array<decltype(ValType() * VecType()), Size>
    operator*(const ValType &val, const std::array<VecType, Size> &a)
{
    std::array<decltype(ValType() * VecType()), Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = a[i] * val;
    return c;
}

template<class ValType, class VecType, size_t Size>
std::array<decltype(VecType() / ValType()), Size>
    operator/(const std::array<VecType, Size> &a, const ValType &val)
{
    std::array<decltype(VecType() / ValType()), Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = a[i] / val;
    return c;
}

template<class ValType, class VecType, size_t Size>
std::array<decltype(VecType() % ValType()), Size>
    operator%(const std::array<VecType, Size> &a, const ValType &val)
{
    std::array<decltype(VecType() % ValType()), Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = a[i] % val;
    return c;
}

template<class LeftType, class RightType, size_t Size>
std::array<decltype(LeftType() & RightType()), Size>
    operator&(const std::array<LeftType, Size> &a, const std::array<LeftType, Size> &b)
{
    std::array<decltype(LeftType() & RightType()), Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = a[i] & b[i];
    return c;
}

template<class LeftType, class RightType, size_t Size>
std::array<decltype(LeftType() | RightType()), Size>
    operator&(const std::array<LeftType, Size> &a, const std::array<LeftType, Size> &b)
{
    std::array<decltype(LeftType() | RightType()), Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = a[i] | b[i];
    return c;
}

template<class LeftType, class RightType, size_t Size>
std::array<decltype(LeftType() ^ RightType()), Size>
    operator&(const std::array<LeftType, Size> &a, const std::array<LeftType, Size> &b)
{
    std::array<decltype(LeftType() ^ RightType()), Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = a[i] ^ b[i];
    return c;
}

template<class Type, size_t Size>
std::array<Type, Size> operator+(const std::array<Type, Size> &a)
{
    std::array<Type, Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = +a[i];
    return c;
}

template<class Type, size_t Size>
std::array<Type, Size> operator-(const std::array<Type, Size> &a)
{
    std::array<Type, Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = -a[i];
    return c;
}

template<class Type, size_t Size>
std::array<Type, Size> operator~(const std::array<Type, Size> &a)
{
    std::array<Type, Size> c;
    for(size_t i = 0; i < Size; ++i)
        c[i] = ~a[i];
    return c;
}

template<class Type, size_t Size>
std::array<Type, Size>& operator++(std::array<Type, Size> &a)
{
    for(size_t i = 0; i < Size; ++i)
        ++a[i];
    return a;
}

template<class Type, size_t Size>
std::array<Type, Size> operator++(std::array<Type, Size> &a, int)
{
    std::array<Type, Size> temp = a;
    for(size_t i = 0; i < Size; ++i)
        a[i]++;
    return temp;
}

template<class Type, size_t Size>
std::array<Type, Size>& operator--(std::array<Type, Size> &a)
{
    for(size_t i = 0; i < Size; ++i)
        --a[i];
    return a;
}

template<class Type, size_t Size>
std::array<Type, Size> operator--(std::array<Type, Size> &a, int)
{
    std::array<Type, Size> temp = a;
    for(size_t i = 0; i < Size; ++i)
        a[i]++;
    return temp;
}

template<class TypeLeft, class TypeRight, size_t Size>
std::array<TypeLeft, Size>& operator+=(std::array<TypeLeft, Size> &a, const std::array<TypeRight, Size> &b)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] += b[i];
    return a;
}

template<class TypeLeft, class TypeRight, size_t Size>
std::array<TypeLeft, Size>& operator-=(std::array<TypeLeft, Size> &a, const std::array<TypeRight, Size> &b)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] -= b[i];
    return a;
}

template<class VecType, class ValType, size_t Size>
std::array<VecType, Size>& operator*=(std::array<VecType, Size> &a, const ValType &val)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] *= val;
    return a;
}

template<class ValType, class VecType, size_t Size>
std::array<VecType, Size>& operator/=(std::array<VecType, Size> &a, const ValType &val)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] /= val;
    return a;
}

template<class ValType, class VecType, size_t Size>
std::array<VecType, Size>& operator%=(std::array<VecType, Size> &a, const ValType &val)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] %= val;
    return a;
}

template<class LeftType, class RightType, size_t Size>
std::array<LeftType, Size>& operator&=(std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] &= b[i];
    return a;
}

template<class LeftType, class RightType, size_t Size>
std::array<LeftType, Size>& operator|=(std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] |= b[i];
    return a;
}

template<class LeftType, class RightType, size_t Size>
std::array<LeftType, Size>& operator^=(std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    for(size_t i = 0; i < Size; ++i)
        a[i] ^= b[i];
    return a;
}

template<class LeftType, class RightType, size_t Size>
inline bool operator==(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    return std::equal(a.cbegin(), a.cend(), b.cbegin());
}

template<class LeftType, class RightType, size_t Size>
inline bool operator!=(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    return !(a == b);
}

template<class LeftType, class RightType, size_t Size>
inline bool operator<(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    return std::lexicographical_compare(a.cbegin(), a.cend(), b.cbegin(), b.cend());
}

template<class LeftType, class RightType, size_t Size>
inline bool operator>(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    return b < a;
}

template<class LeftType, class RightType, size_t Size>
inline bool operator<=(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    return !(a > b);
}

template<class LeftType, class RightType, size_t Size>
inline bool operator>=(const std::array<LeftType, Size> &a, const std::array<RightType, Size> &b)
{
    return !(a < b);
}

#endif