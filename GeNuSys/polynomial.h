#ifndef GENUSYS_ALGEBRA_POLYNOMIAL_H_
#define GENUSYS_ALGEBRA_POLYNOMIAL_H_

namespace GeNuSys
{
    namespace Algebra
    {

        template <typename ElementType>
        class Polynomial
        {
            private:

                unsigned int degree;

                unsigned int capacity;

                ElementType* coefs;

            public:

                Polynomial(unsigned int maxDegree);

                Polynomial(const Polynomial<ElementType>& poly);

                virtual ~Polynomial();

                Polynomial<ElementType>& operator =(const Polynomial<ElementType>& poly);

                int getDegree() const;

                const ElementType& operator[](unsigned int deg) const;
                ElementType& operator[](unsigned int deg);

                void set(unsigned int deg, const ElementType& val);

                ElementType operator()(const ElementType& val) const;

                Polynomial<ElementType> operator +(const Polynomial<ElementType>& poly) const;

                Polynomial<ElementType> operator -(const Polynomial<ElementType>& poly) const;

                Polynomial<ElementType> operator *(const Polynomial<ElementType>& poly) const;

                bool operator ==(const Polynomial<ElementType>& poly);

                bool operator !=(const Polynomial<ElementType>& poly);

        };

    }
}

// Include implementation
#include "polynomial.hpp"

#endif // GENUSYS_ALGEBRA_POLYNOMIAL_H_
