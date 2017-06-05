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

#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys
{

    template<>
    struct ElementTypeTraits<mpf_class>
    {

        typedef mpf_class RealType;

        typedef mpf_class RationalType;

        typedef std::complex<mpf_class> ComplexType;

        typedef mpf_class AbsType;

        typedef mpf_class AbsSqrType;

    };

    template<>
    inline
    const mpf_class& ElementTraits<mpf_class>::zero()
    {
        static const mpf_class zero(0);
        return zero;
    }

    template<>
    inline
    const mpf_class& ElementTraits<mpf_class>::one()
    {
        static const mpf_class one(1);
        return one;
    }

    template<>
    inline
    const mpf_class& ElementTraits<mpf_class>::epsilon()
    {
        static const mpf_class epsilon(ElementTraits<double>::epsilon());
        return epsilon;
    }

    template<>
    template<>
    inline
    mpf_class ElementTraits<mpf_class>::asType<mpf_class>(const mpf_class& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    std::complex<mpf_class> ElementTraits<mpf_class>::asType<std::complex<mpf_class>>(const mpf_class& value)
    {
        return std::complex<mpf_class>(value, mpf_class(0));
    }

    template<>
    template<>
    inline
    mpf_class ElementTraits<int>::asType<mpf_class>(const int& value)
    {
        return mpf_class(value);
    }

    template<>
    template<>
    inline
    std::complex<mpf_class> ElementTraits<int>::asType<std::complex<mpf_class>>(const int& value)
    {
        return std::complex<mpf_class>(mpf_class(value), mpf_class(0));
    }

    template<>
    inline
    mpf_class ElementTraits<mpf_class>::abs(const mpf_class& value)
    {
        mpf_class result;
        mpf_abs(result.get_mpf_t(), value.get_mpf_t());

        return result;
    }

    template<>
    inline
    mpf_class ElementTraits<mpf_class>::absSqr(const mpf_class& value)
    {
        return value * value;
    }

    template<>
    inline
    mpf_class ElementTraits<mpf_class>::pow(const mpf_class& value, unsigned int n)
    {
        mpf_class result;
        mpf_pow_ui(result.get_mpf_t(), value.get_mpf_t(), n);

        return result;
    }

    template<>
    inline
    mpf_class ElementTraits<mpf_class>::sqrt(const mpf_class& value)
    {
        mpf_class result;
        mpf_sqrt(result.get_mpf_t(), value.get_mpf_t());

        return result;
    }

    template<>
    inline
    mpf_class ElementTraits<mpf_class>::div(const mpf_class& a, const mpf_class& b)
    {
        mpf_class result;
        mpf_div(result.get_mpf_t(), a.get_mpf_t(), b.get_mpf_t());

        return result;
    }

}
