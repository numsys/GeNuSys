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
#include <stdexcept>

#include "utils.h"

namespace GeNuSys
{

    template<typename ElementType>
    ExtendedGCD<ElementType>::ExtendedGCD(
        const ElementType& a,
        const ElementType& cA,
        const ElementType& b,
        const ElementType& cB,
        const ElementType& gcd): a(a), cA(cA), b(b), cB(cB), gcd(gcd)
    {
    }

    template<typename ElementType>
    template<typename ConvertedType>
    ConvertedType ElementTraits<ElementType>::asTypeUnsafe(const ElementType& value)
    {
        return ElementTraits<ElementType>::asType<ConvertedType>(value);
    }

    template<typename ElementType>
    ElementType ElementTraits<ElementType>::sgn(const ElementType& value)
    {
        ElementType sign = ElementTraits<ElementType>::zero();
        if (value > ElementTraits<ElementType>::zero())
        {
            sign = ElementTraits<ElementType>::one();
        }
        else if (value < ElementTraits<ElementType>::zero())
        {
            sign = -ElementTraits<ElementType>::one();
        }

        return sign;
    }

    template<typename ElementType>
    ElementType ElementTraits<ElementType>::conj(const ElementType& value)
    {
        return value;
    }

    template<typename ElementType>
    ElementType ElementTraits<ElementType>::pow(const ElementType& value, unsigned int n)
    {
        ElementType result = ElementTraits<ElementType>::one();
        for (ElementType acc = value; n > 0; n /= 2, acc *= acc)
        {
            if (n % 2 == 1)
            {
                result *= acc;
            }
        }

        return result;
    }

    template<typename ElementType>
    ExtendedGCD<ElementType> ElementTraits<ElementType>::egcd(const ElementType& a, const ElementType& b)
    {
        ElementType s = 0, old_s = 1;
        ElementType t = 1, old_t = 0;
        ElementType r = a, old_r = b;
        ElementType tmp, q;
        while (r != 0)
        {
            q = old_r / r;

            tmp = r;
            r = old_r - q * r;
            old_r = tmp;

            tmp = s;
            s = old_s - q * s;
            old_s = tmp;

            tmp = t;
            t = old_t - q * t;
            old_t = tmp;
        }

        ASSERT_EXCEPTION(a * old_t + b * old_s == old_r, std::logic_error);

        return ExtendedGCD<ElementType>(a, old_t, b, old_s, old_r);
    }

}
