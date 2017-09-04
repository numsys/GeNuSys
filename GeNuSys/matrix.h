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

#ifndef GENUSYS_LINALG_MATRIX_H_
#define GENUSYS_LINALG_MATRIX_H_

#include <iostream>
#include <vector>

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        class Vector;

        template<typename ElementType>
        class SparseVector;

        template<typename ElementType>
        class SparseMatrix;

        template<typename ElementType>
        class Matrix
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

                static Matrix<ElementType> diag(const Vector<ElementType>& vct);

                static Matrix<ElementType> diag(const SparseVector<ElementType>& vct);

                static Matrix<ElementType> identity(unsigned int rows, unsigned int cols);

            private:

                unsigned int rows;

                unsigned int cols;

                unsigned int size() const;

                std::vector<ElementType> elem;

                Matrix(unsigned int rows, unsigned int cols, int);

            public:

                Matrix();

                Matrix(unsigned int rows, unsigned int cols);

                template<typename SourceType>
                Matrix(unsigned int rows, unsigned int cols, const std::vector<SourceType>& data);

                Matrix(const Matrix<ElementType>& mat);

                template<typename SourceType>
                Matrix(const Matrix<SourceType>& mat);

                template<typename SourceType>
                Matrix(const SparseMatrix<SourceType>& mat);

                virtual ~Matrix();

                Matrix<ElementType>& operator =(const Matrix<ElementType>& mat);

                template<typename SourceType>
                Matrix<ElementType>& operator =(const Matrix<SourceType>& mat);

                template<typename SourceType>
                Matrix<ElementType>& operator =(const SparseMatrix<SourceType>& mat);

                // Component access

                unsigned int getRows() const;

                unsigned int getCols() const;

                void set(unsigned int i, unsigned int j, const ElementType& value);

                const ElementType& operator()(unsigned int i, unsigned int j) const;

                // Matrix

                Matrix<ElementType> transpose() const;

                Matrix<ElementType> conjugateTranspose() const;

                // Matrix - Value

                Matrix<ElementType> operator *(const ElementType& value) const;

                Matrix<typename ElementTraits<ElementType>::RationalType> operator /(const ElementType& value) const;

                // Matrix - Vector

                Vector<ElementType> operator *(const Vector<ElementType>& vct) const;

                SparseVector<ElementType> operator *=(const Vector<ElementType>& vct) const;

                // Matrix - Sparse Vector

                Vector<ElementType> operator *(const SparseVector<ElementType>& vct) const;

                SparseVector<ElementType> operator *=(const SparseVector<ElementType>& vct) const;

                // Matrix - Matrix

                Matrix<ElementType> operator +(const Matrix<ElementType>& mat) const;

                Matrix<ElementType> operator -(const Matrix<ElementType>& mat) const;

                Matrix<ElementType> operator *(const Matrix<ElementType>& mat) const;

                SparseMatrix<ElementType> operator *=(const Matrix<ElementType>& mat) const;

                bool operator ==(const Matrix<ElementType>& mat);

                bool operator !=(const Matrix<ElementType>& mat);

                // Matrix - Sparse Matrix

                Matrix<ElementType> operator +(const SparseMatrix<ElementType>& mat) const;

                Matrix<ElementType> operator -(const SparseMatrix<ElementType>& mat) const;

                Matrix<ElementType> operator *(const SparseMatrix<ElementType>& mat) const;

                SparseMatrix<ElementType> operator *=(const SparseMatrix<ElementType>& mat) const;

                bool operator ==(const SparseMatrix<ElementType>& mat);

                bool operator !=(const SparseMatrix<ElementType>& mat);

                // IO

                friend std::ostream& operator <<(std::ostream& os, const Matrix<ElementType>& mat)
                {
                    os << "[" << std::endl;
                    for (unsigned int i = 0; i < mat.rows; ++i)
                    {
                        os << "  [ ";
                        for (unsigned int j = 0; j < mat.cols; ++j)
                        {
                            os << " " << mat(i, j);
                        }
                        os << " ]" << std::endl;
                    }
                    os << "]";
                    return os;
                }

        };

    }
}

// Include implementation
#include "matrix.hpp"

#endif // GENUSYS_LINALG_MATRIX_H_
