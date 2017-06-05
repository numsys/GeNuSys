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

#ifndef GENUSYS_LINALG_OPERATIONS_H_
#define GENUSYS_LINALG_OPERATIONS_H_

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        struct Operations
        {

            template<typename ElementType>
            static void vct_mul(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_mul(Vector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_mul(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void vct_mul(SparseVector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_idiv(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_idiv(Vector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_idiv(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void vct_idiv(SparseVector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_mod(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_mod(Vector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_mod(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void vct_mod(SparseVector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_mods(const Vector<ElementType>& vct, const ElementType& value, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_mods(Vector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_mods(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void vct_mods(SparseVector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_div(const Vector<ElementType>& vct, const ElementType& value, Vector<typename ElementTraits<ElementType>::RationalType>& result);

            template<typename ElementType>
            static void vct_div(Vector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_div(const SparseVector<ElementType>& vct, const ElementType& value, SparseVector<typename ElementTraits<ElementType>::RationalType>& result);

            template<typename ElementType>
            static void vct_div(SparseVector<ElementType>& vct, const ElementType& value);

            template<typename ElementType>
            static void vct_add(const Vector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_add(Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static void vct_add(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_add(Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static void vct_add(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_add(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void vct_sub(const Vector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_sub(Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static void vct_sub(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_sub(Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static void vct_sub(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2, Vector<ElementType>& result);

            template<typename ElementType>
            static void vct_sub(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2, SparseVector<ElementType>& result);

            template<typename ElementType>
            static ElementType vct_mul(const Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static ElementType vct_mul(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static ElementType vct_mul(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static ElementType vct_mul(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_eq(const Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_eq(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_eq(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_eq(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_gt(const Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_gt(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_gt(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_gt(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_geq(const Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_geq(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_geq(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_geq(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_lt(const Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_lt(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_lt(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_lt(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_leq(const Vector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_leq(const Vector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_leq(const SparseVector<ElementType>& op1, const Vector<ElementType>& op2);

            template<typename ElementType>
            static bool vct_leq(const SparseVector<ElementType>& op1, const SparseVector<ElementType>& op2);

            template<typename ElementType>
            static void mat_transpose(const Matrix<ElementType>& mat, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_transpose(const SparseMatrix<ElementType>& mat, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_conjugate_transpose(const Matrix<ElementType>& mat, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_conjugate_transpose(const SparseMatrix<ElementType>& mat, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(Matrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(SparseMatrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_idiv(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_idiv(Matrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_idiv(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_idiv(SparseMatrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_mod(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mod(Matrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_mod(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mod(SparseMatrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_mods(const Matrix<ElementType>& mat, const ElementType& value, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mods(Matrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_mods(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mods(SparseMatrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_div(const Matrix<ElementType>& mat, const ElementType& value, Matrix<typename ElementTraits<ElementType>::RationalType>& result);

            template<typename ElementType>
            static void mat_div(Matrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_div(const SparseMatrix<ElementType>& mat, const ElementType& value, SparseMatrix<typename ElementTraits<ElementType>::RationalType>& result);

            template<typename ElementType>
            static void mat_div(SparseMatrix<ElementType>& mat, const ElementType& value);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& mat, const Vector<ElementType>& vct, Vector<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& mat, const Vector<ElementType>& vct, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& mat, const SparseVector<ElementType>& vct, Vector<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& mat, const SparseVector<ElementType>& vct, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& mat, const Vector<ElementType>& vct, Vector<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& mat, const Vector<ElementType>& vct, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& mat, const SparseVector<ElementType>& vct, Vector<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& mat, const SparseVector<ElementType>& vct, SparseVector<ElementType>& result);

            template<typename ElementType>
            static void mat_add(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_add(Matrix<ElementType>& op1, const Matrix<ElementType>& op2);

            template<typename ElementType>
            static void mat_add(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_add(Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);

            template<typename ElementType>
            static void mat_add(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_add(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_sub(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_sub(Matrix<ElementType>& op1, const Matrix<ElementType>& op2);

            template<typename ElementType>
            static void mat_sub(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_sub(Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);

            template<typename ElementType>
            static void mat_sub(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_sub(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, Matrix<ElementType>& result);

            template<typename ElementType>
            static void mat_mul(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2, SparseMatrix<ElementType>& result);

            template<typename ElementType>
            static bool mat_eq(const Matrix<ElementType>& op1, const Matrix<ElementType>& op2);

            template<typename ElementType>
            static bool mat_eq(const Matrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);

            template<typename ElementType>
            static bool mat_eq(const SparseMatrix<ElementType>& op1, const Matrix<ElementType>& op2);

            template<typename ElementType>
            static bool mat_eq(const SparseMatrix<ElementType>& op1, const SparseMatrix<ElementType>& op2);

        };

    }
}

#include "linalg_operations.hpp"

#endif // GENUSYS_LINALG_OPERATIONS_H_
