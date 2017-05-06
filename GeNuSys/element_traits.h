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
