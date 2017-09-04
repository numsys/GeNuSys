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

namespace GeNuSys
{

    template<typename ElementType>
    struct ElementTypeTraits<std::complex<ElementType>>
    {

        typedef std::complex<typename ElementTypeTraits<ElementType>::RealType> RealType;

        typedef std::complex<typename ElementTypeTraits<ElementType>::RationalType> RationalType;

        typedef std::complex<ElementType> ComplexType;

        typedef typename ElementTypeTraits<ElementType>::RealType AbsType;

        typedef ElementType AbsSqrType;

    };

    template<typename ElementType>
    struct ElementTraits<std::complex<ElementType>> : ElementTypeTraits<std::complex<ElementType>>
    {

        static const std::complex<ElementType>& zero();

        static const std::complex<ElementType>& one();

        template<typename ConvertedType>
        static std::complex<ConvertedType> asType(const std::complex<ElementType>& value)
        {
            return std::complex<ConvertedType>(ElementTraits<ElementType>::template asType<ConvertedType>(value.real()), ElementTraits<ElementType>::template asType<ConvertedType>(value.imag()));
        }

        static std::complex<typename ElementTypeTraits<ElementType>::RealType> sgn(const std::complex<ElementType>& value);

        static typename ElementTypeTraits<std::complex<ElementType>>::AbsType abs(const std::complex<ElementType>& value);

        static typename ElementTypeTraits<std::complex<ElementType>>::AbsSqrType absSqr(const std::complex<ElementType>& value);

        static std::complex<ElementType> conj(const std::complex<ElementType>& value);

        static std::complex<ElementType> pow(const std::complex<ElementType>& value, unsigned int n);

        static typename ElementTypeTraits<std::complex<ElementType>>::RealType sqrt(const std::complex<ElementType>& value);

        static typename ElementTypeTraits<std::complex<ElementType>>::RationalType div(const std::complex<ElementType>& a, const std::complex<ElementType>& b);

    };

    template<typename ElementType>
    inline
    const std::complex<ElementType>& ElementTraits<std::complex<ElementType>>::zero()
    {
        static const std::complex<ElementType> zero(ElementTraits<ElementType>::zero(), ElementTraits<ElementType>::zero());
        return zero;
    }

    template<typename ElementType>
    inline
    const std::complex<ElementType>& ElementTraits<std::complex<ElementType>>::one()
    {
        static const std::complex<ElementType> one(ElementTraits<ElementType>::one(), ElementTraits<ElementType>::zero());
        return one;
    }

    template<typename ElementType>
    std::complex<typename ElementTypeTraits<ElementType>::RealType> ElementTraits<std::complex<ElementType>>::sgn(const std::complex<ElementType>& value)
    {
        typedef typename ElementTypeTraits<ElementType>::RealType RealType;

        RealType l = abs(value);
        if (l == ElementTraits<RealType>::zero())
        {
            return zero();
        }

        return std::complex<RealType>(ElementTraits<ElementType>::div(value.real(), l), ElementTraits<ElementType>::div(value.imag(), l));
    }

    template<typename ElementType>
    typename ElementTypeTraits<std::complex<ElementType>>::AbsType ElementTraits<std::complex<ElementType>>::abs(const std::complex<ElementType>& value)
    {
        return ElementTraits<ElementType>::sqrt(value.real() * value.real() + value.imag() * value.imag());
    }

    template<typename ElementType>
    typename ElementTypeTraits<std::complex<ElementType>>::AbsSqrType ElementTraits<std::complex<ElementType>>::absSqr(const std::complex<ElementType>& value)
    {
        return value.real() * value.real() + value.imag() * value.imag();
    }

    template<typename ElementType>
    std::complex<ElementType> ElementTraits<std::complex<ElementType>>::conj(const std::complex<ElementType>& value)
    {
        return std::conj(value);
    }

    template<typename ElementType>
    typename ElementTypeTraits<std::complex<ElementType>>::RealType ElementTraits<std::complex<ElementType>>::sqrt(const std::complex<ElementType>& value)
    {
        typedef typename ElementTypeTraits<ElementType>::RealType RealType;

        if (ElementTraits<ElementType>::abs(value.imag()) <= ElementTraits<typename ElementTraits<ElementType>::AbsType>::epsilon())
        {
            if (value.real() < 0)
            {
                return std::complex<RealType>(ElementTraits<ElementType>::zero(), ElementTraits<RealType>::sqrt(-value.real()));
            }
            return std::complex<RealType>(ElementTraits<RealType>::sqrt(value.real()), ElementTraits<ElementType>::zero());
        }

        RealType l = abs(value);
        RealType re = ElementTraits<RealType>::sqrt((value.real() + l) / 2);
        RealType im = ElementTraits<RealType>::sgn(value.imag()) * ElementTraits<RealType>::sqrt((-value.real() + l) / 2);
        return std::complex<RealType>(re, im);
    }

    template<typename ElementType>
    typename ElementTypeTraits<std::complex<ElementType>>::RationalType ElementTraits<std::complex<ElementType>>::div(const std::complex<ElementType>& a, const std::complex<ElementType>& b)
    {
        ElementType lenSqr = b.real() * b.real() + b.imag() * b.imag();
        std::complex<ElementType> mul = a * std::conj(b);

        return  std::complex<typename ElementTypeTraits<ElementType>::RationalType>(ElementTraits<ElementType>::div(mul.real(), lenSqr), ElementTraits<ElementType>::div(mul.imag(), lenSqr));
    }

    template <class ElementType>
    bool complex_comparator(const std::complex<ElementType>& a, const std::complex<ElementType>& b)
    {
        return real(a) == real(b) ? imag(a) < imag(b) : real(a) < real(b);
    }
}
