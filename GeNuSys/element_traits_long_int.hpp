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
    struct ElementTypeTraits<long int>
    {

        typedef double RealType;

        typedef double RationalType;

        typedef std::complex<long int> ComplexType;

        typedef long int AbsType;

        typedef long int AbsSqrType;

    };

    template<>
    inline
    const long int& ElementTraits<long int>::zero()
    {
        static const long int zero = 0;
        return zero;
    }

    template<>
    inline
    const long int& ElementTraits<long int>::one()
    {
        static const long int one = 1;
        return one;
    }

    template<>
    template<>
    inline
    long int ElementTraits<long int>::asType<long int>(const long int& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    double ElementTraits<long int>::asType<double>(const long int& value)
    {
        return (double) value;
    }

    template<>
    template<>
    inline
    std::complex<long int> ElementTraits<long int>::asType<std::complex<long int>>(const long int& value)
    {
        return std::complex<long int>(value, 0);
    }

    template<>
    template<>
    inline
    std::complex<double> ElementTraits<long int>::asType<std::complex<double>>(const long int& value)
    {
        return std::complex<double>((double) value, 0);
    }

    template<>
    template<>
    inline
    int ElementTraits<long int>::asTypeUnsafe<int>(const long int& value)
    {
        return (int) value;
    }

    template<>
    inline
    long int ElementTraits<long int>::abs(const long int& value)
    {
        return std::abs(value);
    }

    template<>
    inline
    long int ElementTraits<long int>::absSqr(const long int& value)
    {
        return value * value;
    }

    template<>
    inline
    double ElementTraits<long int>::sqrt(const long int& value)
    {
        return std::sqrt((double) value);
    }

    template<>
    inline
    double ElementTraits<long int>::root(const long int& value, int n)
    {
        return std::pow(value, 1.0 / n);
    }

    template<>
    inline
    double ElementTraits<long int>::div(const long int& a, const long int& b)
    {
        return (double) a / b;
    }

    template<>
    inline
    bool ElementTraits<long int>::divisible(const long int& a, const long int& b)
    {
        return a % b == 0;
    }

    template<>
    inline
    long int ElementTraits<long int>::idiv(const long int& a, const long int& b)
    {
        return a / b;
    }

    template<>
    inline
    long int ElementTraits<long int>::mod(const long int& a, const long int& b)
    {
        return (a % b + b) % b;
    }

    template<>
    inline
    long int ElementTraits<long int>::mods(const long int& a, const long int& b)
    {
        long int m = mod(a, b);
        return (m > b / 2) ? m - b : m;
    }

}
