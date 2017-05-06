#include <stdexcept>

#include "utils.h"

namespace GeNuSys
{
    namespace Algebra
    {

        template <typename ElementType>
        Polynomial<ElementType>::Polynomial(unsigned int maxDegree): degree(0), capacity(maxDegree + 1)
        {
            coefs = new ElementType[capacity];
            for (unsigned int i = 0; i < capacity; ++i)
            {
                coefs[i] = ElementTraits<ElementType>::zero();
            }
        }

        template <typename ElementType>
        Polynomial<ElementType>::Polynomial(const Polynomial<ElementType>& poly): degree(poly.degree), capacity(poly.capacity)
        {
            coefs = new ElementType[capacity];
            for (unsigned int i = 0; i <= degree; ++i)
            {
                coefs[i] = poly.coefs[i];
            }
            for (unsigned int i = degree + 1; i < capacity; ++i)
            {
                coefs[i] = ElementTraits<ElementType>::zero();
            }
        }

        template <typename ElementType>
        Polynomial<ElementType>::~Polynomial()
        {
            delete[] coefs;
        }

        template <typename ElementType>
        Polynomial<ElementType>& Polynomial<ElementType>::operator =(const Polynomial<ElementType>& poly)
        {
            if (this != &poly)
            {
                delete[] coefs;

                degree = poly.degree;
                capacity = poly.capacity;
                coefs = new ElementType[capacity];
                for (unsigned int i = 0; i <= degree; ++i)
                {
                    coefs[i] = poly.coefs[i];
                }
                for (unsigned int i = degree + 1; i < capacity; ++i)
                {
                    coefs[i] = ElementTraits<ElementType>::zero();
                }
            }
            return *this;
        }

        template <typename ElementType>
        int Polynomial<ElementType>::getDegree() const
        {
            return degree;
        }

        template <typename ElementType>
        const ElementType& Polynomial<ElementType>::operator[](unsigned int deg) const
        {
            ASSERT_EXCEPTION(deg < capacity, std::out_of_range);
            return coefs[deg];
        }

        template <typename ElementType>
        ElementType& Polynomial<ElementType>::operator[](unsigned int deg)
        {
            ASSERT_EXCEPTION(deg < capacity, std::out_of_range);
            return coefs[deg];
        }

        template <typename ElementType>
        void Polynomial<ElementType>::set(unsigned int deg, const ElementType& val)
        {
            ASSERT_EXCEPTION(deg < capacity, std::out_of_range);
            coefs[deg] = val;
            if (deg > degree && val != ElementTraits<ElementType>::zero())
            {
                degree = deg;
            }
            else
            {
                while (coefs[degree] == ElementTraits<ElementType>::zero() && degree > 0)
                {
                    --degree;
                }
            }
        }

        template <typename ElementType>
        ElementType Polynomial<ElementType>::operator()(const ElementType& val) const
        {
            ElementType prod = coefs[degree];
            for (int i = degree - 1; i >= 0; --i)
            {
                prod *= val;
                prod += coefs[i];
            }
            return prod;
        }

        template <typename ElementType>
        Polynomial<ElementType> Polynomial<ElementType>::operator +(const Polynomial<ElementType>& poly) const
        {
            if (degree == poly.degree)
            {
                Polynomial<ElementType> result(degree);
                for (unsigned int i = 0; i <= degree; ++i)
                {
                    result[i] = coefs[i] + poly.coefs[i];
                }
                while (coefs[degree] == ElementTraits<ElementType>::zero() && degree > 0)
                {
                    --result.degree;
                }
                return result;
            }
            else if (degree < poly.degree)
            {
                Polynomial<ElementType> result(poly.degree);
                for (unsigned int i = 0; i <= degree; ++i)
                {
                    result.coefs[i] = coefs[i] + poly.coefs[i];
                }
                for (unsigned int i = degree + 1; i <= poly.degree; ++i)
                {
                    result.coefs[i] = poly.coefs[i];
                }
                return result;
            }
            else
            {
                Polynomial<ElementType> result(degree);
                for (unsigned int i = 0; i <= poly.degree; ++i)
                {
                    result.coefs[i] = coefs[i] + poly.coefs[i];
                }
                for (unsigned int i = poly.degree + 1; i <= degree; ++i)
                {
                    result.coefs[i] = coefs[i];
                }
                return result;
            }
        }

        template <typename ElementType>
        Polynomial<ElementType> Polynomial<ElementType>::operator -(const Polynomial<ElementType>& poly) const
        {
            if (degree == poly.degree)
            {
                Polynomial<ElementType> result(degree);
                for (unsigned int i = 0; i <= degree; ++i)
                {
                    result[i] = coefs[i] - poly.coefs[i];
                }
                while (coefs[degree] == ElementTraits<ElementType>::zero() && degree > 0)
                {
                    --result.degree;
                }
                return result;
            }
            else if (degree < poly.degree)
            {
                Polynomial<ElementType> result(poly.degree);
                for (unsigned int i = 0; i <= degree; ++i)
                {
                    result.coefs[i] = coefs[i] - poly.coefs[i];
                }
                for (unsigned int i = degree + 1; i <= poly.degree; ++i)
                {
                    result.coefs[i] = -poly.coefs[i];
                }
                return result;
            }
            else
            {
                Polynomial<ElementType> result(degree);
                for (unsigned int i = 0; i <= poly.degree; ++i)
                {
                    result.coefs[i] = coefs[i] - poly.coefs[i];
                }
                for (unsigned int i = poly.degree + 1; i <= degree; ++i)
                {
                    result.coefs[i] = coefs[i];
                }
                return result;
            }
        }

        template <typename ElementType>
        Polynomial<ElementType> Polynomial<ElementType>::operator *(const Polynomial<ElementType>& poly) const
        {
            Polynomial<ElementType> result(degree + poly.degree);
            for (unsigned int i = 0; i <= result.degree; ++i)
            {
                result.coefs[i] = ElementTraits<ElementType>::zero();
            }
            for (unsigned int i = 0; i <= degree; ++i)
            {
                for (unsigned int j = 0; j <= poly.degree; ++j)
                {
                    result.coefs[i + j] += coefs[i] * poly.coefs[j];
                }
            }
            return result;
        }

        template <typename ElementType>
        bool Polynomial<ElementType>::operator ==(const Polynomial<ElementType>& poly)
        {
            if (degree != poly.degree)
            {
                return false;
            }
            bool l = true;
            for (unsigned int i = 0; i <= degree && l; ++i)
            {
                l = (coefs[i] == poly.coefs[i]);
            }
            return l;
        }

        template <typename ElementType>
        bool Polynomial<ElementType>::operator !=(const Polynomial<ElementType>& poly)
        {
            return !(*this == poly);
        }

    }
}
