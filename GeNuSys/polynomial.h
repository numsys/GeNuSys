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
