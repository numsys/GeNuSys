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
    struct ElementTypeTraits<long long>
    {

        typedef double RealType;

        typedef double RationalType;

        typedef std::complex<long long> ComplexType;

        typedef long long AbsType;

        typedef long long AbsSqrType;

    };

    template<>
    inline
    const long long& ElementTraits<long long>::zero()
    {
        static const long long zero = 0;
        return zero;
    }

    template<>
    inline
    const long long& ElementTraits<long long>::one()
    {
        static const long long one = 1;
        return one;
    }

    template<>
    template<>
    inline
    long long ElementTraits<long long>::asType<long long>(const long long& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    double ElementTraits<long long>::asType<double>(const long long& value)
    {
        return (double) value;
    }

    template<>
    template<>
    inline
    std::complex<long long> ElementTraits<long long>::asType<std::complex<long long>>(const long long& value)
    {
        return std::complex<long long>(value, 0);
    }

    template<>
    template<>
    inline
    std::complex<double> ElementTraits<long long>::asType<std::complex<double>>(const long long& value)
    {
        return std::complex<double>((double) value, 0);
    }

    template<>
    template<>
    inline
    unsigned long ElementTraits<long long>::asTypeUnsafe<unsigned long>(const long long& value)
    {
        return (unsigned long) value;
    }

    template<>
    template<>
    inline
    int ElementTraits<long long>::asTypeUnsafe<>(const long long& value)
    {
        return (int) value;
    }

    template<>
    template<>
    inline
    long ElementTraits<long long>::asTypeUnsafe<long>(const long long& value)
    {
        return (long) value;
    }
    
    template<>
    inline
    long long ElementTraits<long long>::abs(const long long& value)
    {
        return std::abs(value);
    }

    template<>
    inline
    long long ElementTraits<long long>::absSqr(const long long& value)
    {
        return value * value;
    }

    template<>
    inline
    double ElementTraits<long long>::sqrt(const long long& value)
    {
        return std::sqrt((double) value);
    }

    template<>
    inline
    double ElementTraits<long long>::root(const long long& value, int n)
    {
        return std::pow(value, 1.0 / n);
    }

    template<>
    inline
    double ElementTraits<long long>::div(const long long& a, const long long& b)
    {
        return (double) a / b;
    }

    template<>
    inline
    bool ElementTraits<long long>::divisible(const long long& a, const long long& b)
    {
        return a % b == 0;
    }

    template<>
    inline
    long long ElementTraits<long long>::idiv(const long long& a, const long long& b)
    {
        return a / b;
    }

    template<>
    inline
    long long ElementTraits<long long>::mod(const long long& a, const long long& b)
    {
        return (a % b + b) % b;
    }

    template<>
    inline
    long long ElementTraits<long long>::mods(const long long& a, const long long& b)
    {
        long long m = mod(a, b);
        return (m > b / 2) ? m - b : m;
    }

}
