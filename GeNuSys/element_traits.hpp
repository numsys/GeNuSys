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
