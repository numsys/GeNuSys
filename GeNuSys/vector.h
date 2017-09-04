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

#ifndef GENUSYS_LINALG_VECTOR_H_
#define GENUSYS_LINALG_VECTOR_H_

#include <iostream>
#include <vector>

#include "element_traits.h"

#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        class SparseVector;

        template<typename ElementType>
        class Matrix;

        template<typename ElementType>
        class SparseMatrix;

        template<typename ElementType>
        class Vector
        {

                friend struct Traits;

                friend struct Algorithms;

                friend struct Operations;

                template<typename T>
                friend class Vector;

                template<typename T>
                friend class SparseVector;

                template<typename T>
                friend class Matrix;

                template<typename T>
                friend class SparseMatrix;

                template<unsigned int P>
                friend struct PNorm;

                template<typename T>
                friend class OperatorNorm;

            private:

                unsigned int length;

                std::vector<ElementType> elem;

                Vector(unsigned int length, int);

            public:

                Vector();

                Vector(unsigned int length);

                Vector(const Vector<ElementType>& vct);

                template<typename SourceType>
                Vector(const Vector<SourceType>& vct);

                template<typename SourceType>
                Vector(const SparseVector<SourceType>& vct);

                virtual ~Vector();

                Vector<ElementType>& operator =(const Vector<ElementType>& vct);

                template<typename SourceType>
                Vector<ElementType>& operator =(const Vector<SourceType>& vct);

                template<typename SourceType>
                Vector<ElementType>& operator =(const SparseVector<SourceType>& vct);

                // Component access

                unsigned int getLength() const;

                const ElementType& operator [](unsigned int idx) const;

                void set(unsigned int idx, const ElementType& value);

                // Vector - Value

                Vector<ElementType> operator *(const ElementType& value) const;

                Vector<typename ElementTraits<ElementType>::RationalType> operator /(const ElementType& value) const;

                // Vector - Vector

                ElementType operator *(const Vector<ElementType>& vct) const;

                Vector<ElementType> operator +(const Vector<ElementType>& vct) const;

                Vector<ElementType> operator -(const Vector<ElementType>& vct) const;

                bool operator ==(const Vector<ElementType>& vct) const;

                bool operator !=(const Vector<ElementType>& vct) const;

                bool operator <(const Vector<ElementType>& vct) const;

                bool operator <=(const Vector<ElementType>& vct) const;

                bool operator >(const Vector<ElementType>& vct) const;

                bool operator >=(const Vector<ElementType>& vct) const;

                // Vector - Sparse Vector

                ElementType operator *(const SparseVector<ElementType>& vct) const;

                Vector<ElementType> operator +(const SparseVector<ElementType>& vct) const;

                Vector<ElementType> operator -(const SparseVector<ElementType>& vct) const;

                bool operator ==(const SparseVector<ElementType>& vct) const;

                bool operator !=(const SparseVector<ElementType>& vct) const;

                bool operator <(const SparseVector<ElementType>& vct) const;

                bool operator <=(const SparseVector<ElementType>& vct) const;

                bool operator >(const SparseVector<ElementType>& vct) const;

                bool operator >=(const SparseVector<ElementType>& vct) const;

                // IO

                friend std::ostream& operator <<(std::ostream& os, const Vector<ElementType>& vct)
                {
                    os << "[ ";
                    for (unsigned int i = 0; i < vct.length; ++i)
                    {
                        os << vct.elem[i] << " ";
                    }
                    os << "]";
                    return os;
                }

        };

    }
}

// Include implementation
#include "vector.hpp"

#endif // GENUSYS_LINALG_VECTOR_H_
