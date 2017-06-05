/*
GeNuSys - computations with generalized number systems
Copyright (C) 2015-2017  Bence Németh
Copyright (C) 2017  Tamás Krutki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <complex>

namespace GeNuSys
{

    template<>
    struct ElementTypeTraits<int>
    {

        typedef double RealType;

        typedef double RationalType;

        typedef std::complex<int> ComplexType;

        typedef int AbsType;

        typedef int AbsSqrType;

    };

    template<>
    inline
    const int& ElementTraits<int>::zero()
    {
        static const int zero = 0;
        return zero;
    }

    template<>
    inline
    const int& ElementTraits<int>::one()
    {
        static const int one = 1;
        return one;
    }

    template<>
    template<>
    inline
    int ElementTraits<int>::asType<int>(const int& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    long long ElementTraits<int>::asType<long long>(const int& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    long ElementTraits<int>::asType<long>(const int& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    unsigned long ElementTraits<int>::asTypeUnsafe<unsigned long>(const int& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    double ElementTraits<int>::asType<double>(const int& value)
    {
        return (double) value;
    }

    template<>
    template<>
    inline
    std::complex<int> ElementTraits<int>::asType<std::complex<int>>(const int& value)
    {
        return std::complex<int>(value, 0);
    }

    template<>
    template<>
    inline
    std::complex<long long> ElementTraits<int>::asType<std::complex<long long>>(const int& value)
    {
        return std::complex<long long>(value, 0);
    }

    template<>
    template<>
    inline
    std::complex<double> ElementTraits<int>::asType<std::complex<double>>(const int& value)
    {
        return std::complex<double>((double) value, 0);
    }

    template<>
    inline
    int ElementTraits<int>::abs(const int& value)
    {
        return std::abs(value);
    }

    template<>
    inline
    int ElementTraits<int>::absSqr(const int& value)
    {
        return value * value;
    }

    template<>
    inline
    double ElementTraits<int>::sqrt(const int& value)
    {
        return std::sqrt((double) value);
    }

    template<>
    inline
    double ElementTraits<int>::root(const int& value, int n)
    {
        return std::pow(value, 1.0 / n);
    }

    template<>
    inline
    double ElementTraits<int>::div(const int& a, const int& b)
    {
        return (double) a / b;
    }

    template<>
    inline
    bool ElementTraits<int>::divisible(const int& a, const int& b)
    {
        return a % b == 0;
    }

    template<>
    inline
    int ElementTraits<int>::idiv(const int& a, const int& b)
    {
        return a / b;
    }

    template<>
    inline
    int ElementTraits<int>::mod(const int& a, const int& b)
    {
        return (a % b + b) % b;
    }

    template<>
    inline
    int ElementTraits<int>::mods(const int& a, const int& b)
    {
        int m = mod(a, b);
        return (m > b / 2) ?  m - b : m;
    }

}
