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

#ifndef GENUSYS_ELEMENT_TRAITS_H_
#define GENUSYS_ELEMENT_TRAITS_H_

namespace GeNuSys
{

    template<typename ElementType>
    struct ExtendedGCD
    {

        const ElementType a;

        const ElementType cA;

        const ElementType b;

        const ElementType cB;

        const ElementType gcd;

        ExtendedGCD(const ElementType& a, const ElementType& cA, const ElementType& b, const ElementType& cB, const ElementType& gcd);

    };

    template<typename ElementType>
    struct ElementTypeTraits
    {

    };

    template<typename ElementType>
    struct ElementTraits : ElementTypeTraits<ElementType>
    {

        static const ElementType& zero();

        static const ElementType& one();

        static const ElementType& epsilon();

        template<typename ConvertedType>
        static ConvertedType asType(const ElementType& value);

        template<typename ConvertedType>
        static ConvertedType asTypeUnsafe(const ElementType& value);

        static ElementType sgn(const ElementType& value);

        static typename ElementTypeTraits<ElementType>::AbsType abs(const ElementType& value);

        static typename ElementTypeTraits<ElementType>::AbsSqrType absSqr(const ElementType& value);

        static ElementType conj(const ElementType& value);

        static ElementType pow(const ElementType& value, unsigned int n);

        static typename ElementTypeTraits<ElementType>::RealType sqrt(const ElementType& value);

        static typename ElementTypeTraits<ElementType>::RealType root(const ElementType& value, int n);

        static typename ElementTypeTraits<ElementType>::RationalType div(const ElementType& a, const ElementType& b);

        static bool divisible(const ElementType& a, const ElementType& b);

        static ElementType idiv(const ElementType& a, const ElementType& b);

        static ElementType mod(const ElementType& a, const ElementType& b);

        static ElementType mods(const ElementType& a, const ElementType& b);

        static ExtendedGCD<ElementType> egcd(const ElementType& a, const ElementType& b);

    };

}

// Include implementation
#include "element_traits.hpp"
#include "element_traits_int.hpp"
#include "element_traits_long_long.hpp"
#include "element_traits_long_int.hpp"
#include "element_traits_double.hpp"
#include "element_traits_complex.hpp"
#ifdef __unix__
#include "element_traits_mpz.hpp"
#include "element_traits_mpq.hpp"
#include "element_traits_mpf.hpp"
#endif // __unix__

#endif // GENUSYS_ELEMENT_TRAITS_H_
