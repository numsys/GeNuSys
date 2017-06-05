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

#ifndef GENUSYS_LINALG_SPARSE_MATRIX_H_
#define GENUSYS_LINALG_SPARSE_MATRIX_H_

#include <iostream>
#include <vector>

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        class Vector;

        template<typename ElementType>
        class SparseVector;

        template<typename ElementType>
        class Matrix;

        template<typename ElementType>
        class SparseMatrix
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

                template<unsigned int p>
                friend struct PNorm;

                friend struct FrobeniusNorm;

                template<typename T>
                friend class OperatorNorm;

            public:

                static SparseMatrix<ElementType> diag(const Vector<ElementType>& vct);

                static SparseMatrix<ElementType> diag(const SparseVector<ElementType>& vct);

                static SparseMatrix<ElementType> identity(unsigned int rows, unsigned int cols);

            public:

                struct Entry
                {

                    unsigned int col_idx;

                    ElementType value;

                    Entry();

                    Entry(unsigned int col_idx, const ElementType& value);

                };

            private:

                unsigned int rows;

                unsigned int cols;

                std::vector<unsigned int> row_ptr;

                std::vector<Entry> elem;

                unsigned int size() const;

                unsigned int search(unsigned int row_idx, unsigned int col_idx) const;

                void push(unsigned int col_idx, const ElementType& value);

                SparseMatrix(unsigned int rows, unsigned int cols, unsigned int capacity, unsigned int size);

                SparseMatrix(unsigned int rows, unsigned int cols, unsigned int capacity);

            public:

                SparseMatrix();

                SparseMatrix(unsigned int rows, unsigned int cols);

                template<typename SourceType>
                SparseMatrix(unsigned int rows, unsigned int cols, const std::vector<SourceType>& data);

                SparseMatrix(const SparseMatrix<ElementType>& mat);

                template<typename SourceType>
                SparseMatrix(const SparseMatrix<SourceType>& mat);

                template<typename SourceType>
                SparseMatrix(const Matrix<SourceType>& mat);

                virtual ~SparseMatrix();

                SparseMatrix<ElementType>& operator =(const SparseMatrix<ElementType>& mat);

                template<typename SourceType>
                SparseMatrix<ElementType>& operator =(const SparseMatrix<SourceType>& mat);

                template<typename SourceType>
                SparseMatrix<ElementType>& operator =(const Matrix<SourceType>& mat);

                // Component access

                unsigned int getRows() const;

                unsigned int getCols() const;

                void set(unsigned int row_idx, unsigned int col_idx, const ElementType& value);

                const ElementType& operator()(unsigned int row_idx, unsigned int col_idx) const;

                // Matrix

                SparseMatrix<ElementType> transpose() const;

                SparseMatrix<ElementType> conjugateTranspose() const;

                // Sparse Matrix - Value

                SparseMatrix<ElementType> operator *(const ElementType& value) const;

                SparseMatrix<typename ElementTraits<ElementType>::RationalType> operator /(const ElementType& value) const;

                // Sparse Matrix - Vector

                Vector<ElementType> operator *(const Vector<ElementType>& vct) const;

                SparseVector<ElementType> operator *=(const Vector<ElementType>& vct) const;

                // Sparse Matrix - Sparse Vector

                Vector<ElementType> operator *(const SparseVector<ElementType>& vct) const;

                SparseVector<ElementType> operator *=(const SparseVector<ElementType>& vct) const;

                // Sparse Matrix - Matrix

                Matrix<ElementType> operator +(const Matrix<ElementType>& mat) const;

                Matrix<ElementType> operator -(const Matrix<ElementType>& mat) const;

                Matrix<ElementType> operator *(const Matrix<ElementType>& mat) const;

                SparseMatrix<ElementType> operator *=(const Matrix<ElementType>& mat) const;

                bool operator ==(const Matrix<ElementType>& mat);

                bool operator !=(const Matrix<ElementType>& mat);

                // Sparse Matrix - Sparse Matrix

                SparseMatrix<ElementType> operator +(const SparseMatrix<ElementType>& mat) const;

                SparseMatrix<ElementType> operator -(const SparseMatrix<ElementType>& mat) const;

                Matrix<ElementType> operator *(const SparseMatrix<ElementType>& mat) const;

                SparseMatrix<ElementType> operator *=(const SparseMatrix<ElementType>& mat) const;

                bool operator ==(const SparseMatrix<ElementType>& mat);

                bool operator !=(const SparseMatrix<ElementType>& mat);

                // IO

                friend std::ostream& operator <<(std::ostream& os, const SparseMatrix<ElementType>& mat)
                {
                    os << "[ " << std::endl;
                    for (unsigned int row = 0; row < mat.rows; ++row)
                    {
                        for (unsigned int i = mat.row_ptr[row]; i < mat.row_ptr[row + 1]; ++i)
                        {
                            os << " (" << row << ", " << mat.elem[i].col_idx << ") -> " << mat.elem[i].value << std::endl;
                        }
                    }
                    os << "]";
                    return os;
                }

        };

    }
}

// Include implementation
#include "sparse_matrix.hpp"

#endif // GENUSYS_LINALG_SPARSE_MATRIX_H_
