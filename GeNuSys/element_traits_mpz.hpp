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

#include <stdexcept>
#include <complex>

#include <gmp.h>
#include <gmpxx.h>

#include "utils.h"

namespace GeNuSys
{

    template<>
    struct ElementTypeTraits<mpz_class>
    {

        typedef mpf_class RealType;

        typedef mpq_class RationalType;

        typedef std::complex<mpz_class> ComplexType;

        typedef mpz_class AbsType;

        typedef mpz_class AbsSqrType;

    };

    template<>
    inline
    const mpz_class& ElementTraits<mpz_class>::zero()
    {
        static const mpz_class zero(0);
        return zero;
    }

    template<>
    inline
    const mpz_class& ElementTraits<mpz_class>::one()
    {
        static const mpz_class one(1);
        return one;
    }

    template<>
    template<>
    inline
    mpz_class ElementTraits<mpz_class>::asType<mpz_class>(const mpz_class& value)
    {
        return value;
    }

    template<>
    template<>
    inline
    mpf_class ElementTraits<mpz_class>::asType<mpf_class>(const mpz_class& value)
    {
        return mpf_class(value);
    }

    template<>
    template<>
    inline
    mpq_class ElementTraits<mpz_class>::asType<mpq_class>(const mpz_class& value)
    {
        return mpq_class(value);
    }

    template<>
    template<>
    inline
    int ElementTraits<mpz_class>::asTypeUnsafe<int>(const mpz_class& value)
    {
        if (value.fits_sint_p())
        {
            return value.get_si();
        }
        else
        {
            throw std::out_of_range{"Conversion to int failed: out of range"};
        }
    }
    
    template<>
    template<>
    inline
    long int ElementTraits<mpz_class>::asTypeUnsafe<long int>(const mpz_class& value)
    {
        if (value.fits_slong_p())
        {
            return value.get_si();
        }
        else
        {
            throw std::out_of_range{"Conversion to long int failed: out of range"};
        }
    }

    template<>
    template<>
    inline
    long int ElementTraits<mpz_class>::asType<long int>(const mpz_class& value)
    {
        return ElementTraits<mpz_class>::asTypeUnsafe<long int>(value);
    }

    template<>
    template<>
    inline
    long long ElementTraits<mpz_class>::asType<long long>(const mpz_class& value)
    {
        return ElementTraits<mpz_class>::asTypeUnsafe<long int>(value);
    }
    
    template<>
    template<>
    inline
    unsigned long int ElementTraits<mpz_class>::asTypeUnsafe<unsigned long int>(const mpz_class& value)
    {
        if (value.fits_ulong_p())
        {
            return value.get_ui();
        }
        else
        {
            throw std::out_of_range{"Conversion to unsigned long int failed: out of range"};
        }
    }

    template<>
    template<>
    inline
    std::complex<mpz_class> ElementTraits<mpz_class>::asType<std::complex<mpz_class>>(const mpz_class& value)
    {
        return std::complex<mpz_class>(value, mpz_class(0));
    }

    template<>
    template<>
    inline
    std::complex<mpf_class> ElementTraits<mpz_class>::asType<std::complex<mpf_class>>(const mpz_class& value)
    {
        return std::complex<mpf_class>(value, mpz_class(0));
    }

    template<>
    template<>
    inline
    std::complex<mpq_class> ElementTraits<mpz_class>::asType<std::complex<mpq_class>>(const mpz_class& value)
    {
        return std::complex<mpq_class>(value, mpz_class(0));
    }

    template<>
    template<>
    inline
    mpz_class ElementTraits<int>::asType<mpz_class>(const int& value)
    {
        return mpz_class(value);
    }

    template<>
    template<>
    inline
    std::complex<mpz_class> ElementTraits<int>::asType<std::complex<mpz_class>>(const int& value)
    {
        return std::complex<mpz_class>(mpz_class(value), mpz_class(0));
    }

    template<>
    template<>
    inline
    mpz_class ElementTraits<long long>::asType<mpz_class>(const long long& value)
    {
        return mpz_class(std::to_string(value));
    }

    template<>
    template<>
    inline
    mpz_class ElementTraits<long>::asType<mpz_class>(const long& value)
    {
        return mpz_class(value);
    }
    
    template<>
    template<>
    inline
    std::complex<mpz_class> ElementTraits<long long>::asType<std::complex<mpz_class>>(const long long& value)
    {
        return std::complex<mpz_class>(mpz_class(std::to_string(value)), mpz_class(0));
    }

    template<>
    template<>
    inline
    mpz_class ElementTraits<double>::asType<mpz_class>(const double&)
    {
        ASSERT_EXCEPTION(!"NOT IMPLEMENTED", std::runtime_error);
    }

    template<>
    template<>
    inline
    std::complex<mpz_class> ElementTraits<double>::asType<std::complex<mpz_class>>(const double&)
    {
        ASSERT_EXCEPTION(!"NOT IMPLEMENTED", std::runtime_error);
    }

    template<>
    inline
    mpz_class ElementTraits<mpz_class>::abs(const mpz_class& value)
    {
        mpz_class result;
        mpz_abs(result.get_mpz_t(), value.get_mpz_t());

        return result;
    }

    template<>
    inline
    mpz_class ElementTraits<mpz_class>::absSqr(const mpz_class& value)
    {
        return value * value;
    }

    template<>
    inline
    mpz_class ElementTraits<mpz_class>::pow(const mpz_class& value, unsigned int n)
    {
        mpz_class result;
        mpz_pow_ui(result.get_mpz_t(), value.get_mpz_t(), n);

        return result;
    }

    template<>
    inline
    mpf_class ElementTraits<mpz_class>::sqrt(const mpz_class& value)
    {
        mpf_class result;
        mpf_sqrt(result.get_mpf_t(), mpf_class(value).get_mpf_t());

        return result;
    }

    template<>
    inline
    mpq_class ElementTraits<mpz_class>::div(const mpz_class& a, const mpz_class& b)
    {
        mpq_class result;
        mpq_set_num(result.get_mpq_t(), a.get_mpz_t());
        mpq_set_den(result.get_mpq_t(), b.get_mpz_t());
        mpq_canonicalize(result.get_mpq_t());

        return result;
    }

    template<>
    inline
    bool ElementTraits<mpz_class>::divisible(const mpz_class& a, const mpz_class& b)
    {
        return mpz_divisible_p(a.get_mpz_t(), b.get_mpz_t()) != 0;
    }

    template<>
    inline
    mpz_class ElementTraits<mpz_class>::idiv(const mpz_class& a, const mpz_class& b)
    {
        mpz_class result;
        mpz_fdiv_q(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

        return result;
    }

    template<>
    inline
    mpz_class ElementTraits<mpz_class>::mod(const mpz_class& a, const mpz_class& b)
    {
        mpz_class result;
        mpz_mod(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

        return result;
    }

    template<>
    inline
    mpz_class ElementTraits<mpz_class>::mods(const mpz_class& a, const mpz_class& b)
    {
        mpz_class result;
        mpz_mod(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

        return (result > b / 2) ? result - b : result;
    }

    template<>
    inline
    ExtendedGCD<mpz_class> ElementTraits<mpz_class>::egcd(const mpz_class& a, const mpz_class& b)
    {
        mpz_class gcd, aC, bC;
        mpz_gcdext(gcd.get_mpz_t(), aC.get_mpz_t(), bC.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

        return ExtendedGCD<mpz_class>(a, aC, b, bC, gcd);
    }

}
