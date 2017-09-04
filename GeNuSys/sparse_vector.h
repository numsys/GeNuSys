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

#ifndef GENUSYS_LINALG_SPARSE_VECTOR_H_
#define GENUSYS_LINALG_SPARSE_VECTOR_H_

#include <iostream>
#include <vector>

#include "element_traits.h"

#include "vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        class Vector;

        template<typename ElementType>
        class Matrix;

        template<typename ElementType>
        class SparseMatrix;

        template<typename ElementType>
        class SparseVector
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

            public:

                struct Entry
                {

                    unsigned int idx;

                    ElementType value;

                    Entry();

                    Entry(unsigned int idx, const ElementType& value);

                };

            private:

                unsigned int length;

                std::vector<Entry> elem;

                unsigned int size() const;

                unsigned int search(unsigned int idx) const;

                void push(unsigned int idx, const ElementType& value);

                SparseVector(unsigned int length, unsigned int size);

            public:

                SparseVector();

                SparseVector(unsigned int length);

                SparseVector(const SparseVector<ElementType>& vct);

                template<typename SourceType>
                SparseVector(const SparseVector<SourceType>& vct);

                template<typename SourceType>
                SparseVector(const Vector<SourceType>& vct);

                virtual ~SparseVector();

                SparseVector<ElementType>& operator =(const SparseVector<ElementType>& vct);

                template<typename SourceType>
                SparseVector<ElementType>& operator =(const SparseVector<SourceType>& vct);

                template<typename SourceType>
                SparseVector<ElementType>& operator =(const Vector<SourceType>& vct);

                // Component access

                unsigned int getLength() const;

                void set(unsigned int idx, const ElementType& value);

                const ElementType& operator [](unsigned int idx) const;

                // Vector - Value

                SparseVector<ElementType> operator *(const ElementType& value) const;

                SparseVector<typename ElementTraits<ElementType>::RationalType> operator /(const ElementType& value) const;

                // Sparse Vector - Vector

                ElementType operator *(const Vector<ElementType>& vct) const;

                Vector<ElementType> operator +(const Vector<ElementType>& vct) const;

                Vector<ElementType> operator -(const Vector<ElementType>& vct) const;

                bool operator ==(const Vector<ElementType>& vct) const;

                bool operator !=(const Vector<ElementType>& vct) const;

                bool operator <(const Vector<ElementType>& vct) const;

                bool operator <=(const Vector<ElementType>& vct) const;

                bool operator >(const Vector<ElementType>& vct) const;

                bool operator >=(const Vector<ElementType>& vct) const;

                // Sparse Vector - Sparse Vector

                ElementType operator *(const SparseVector<ElementType>& vct) const;

                SparseVector<ElementType> operator +(const SparseVector<ElementType>& vct) const;

                SparseVector<ElementType> operator -(const SparseVector<ElementType>& vct) const;

                bool operator ==(const SparseVector<ElementType>& vct) const;

                bool operator !=(const SparseVector<ElementType>& vct) const;

                bool operator <(const SparseVector<ElementType>& vct) const;

                bool operator <=(const SparseVector<ElementType>& vct) const;

                bool operator >(const SparseVector<ElementType>& vct) const;

                bool operator >=(const SparseVector<ElementType>& vct) const;

                // IO

                friend std::ostream& operator <<(std::ostream& os, const SparseVector<ElementType>& vct)
                {
                    os << "[";
                    for (unsigned int i = 0; i < vct.size(); ++i)
                    {
                        os << " (" << vct.elem[i].idx << ", " << vct.elem[i].value << ")";
                    }
                    os << " ]";
                    return os;
                }

        };

    }
}

// Include implementation
#include "sparse_vector.hpp"

#endif // GENUSYS_LINALG_SPARSE_VECTOR_H_
