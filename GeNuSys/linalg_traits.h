#ifndef GENUSYS_LINALG_LINALG_TRAITS_H_
#define GENUSYS_LINALG_LINALG_TRAITS_H_

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        struct Traits
        {

            template<typename SourceType, typename ConvertedType>
            static Matrix<ConvertedType> convertUnsafe(const Matrix<SourceType>& mat);

            template<typename SourceType, typename ConvertedType>
            static SparseMatrix<ConvertedType> convertUnsafe(const SparseMatrix<SourceType>& mat);

            template<typename SourceType, typename ConvertedType>
            static Vector<ConvertedType> convertUnsafe(const Vector<SourceType>& vct);

            template<typename SourceType, typename ConvertedType>
            static SparseVector<ConvertedType> convertUnsafe(const SparseVector<SourceType>& vct);

            template<typename ElementType>
            static void swapRows(Matrix<ElementType>& mat, const int rowA, const int rowB);

            template<typename ElementType>
            static void swapCols(Matrix<ElementType>& mat, const int colA, const int colB);

            template<typename ElementType>
            static Vector<ElementType> getSubVector(const Vector<ElementType>& vct, unsigned int from, unsigned int to);

            template<typename ElementType>
            static SparseVector<ElementType> getSubVector(const SparseVector<ElementType>& vct, unsigned int from, unsigned int to);

            template<typename ElementType>
            static Vector<ElementType> getRow(const Matrix<ElementType>& mat, unsigned int row);

            template<typename ElementType>
            static SparseVector<ElementType> getRow(const SparseMatrix<ElementType>& mat, unsigned int row);

            template<typename ElementType>
            static Vector<ElementType> getRow(const Matrix<ElementType>& mat, unsigned int row, unsigned int from, unsigned int to);

            template<typename ElementType>
            static SparseVector<ElementType> getRow(const SparseMatrix<ElementType>& mat, unsigned int row, unsigned int from, unsigned int to);

            template<typename ElementType>
            static Matrix<ElementType> getRows(const Matrix<ElementType>& mat, unsigned int from, unsigned int to);

            template<typename ElementType>
            static SparseMatrix<ElementType> getRows(const SparseMatrix<ElementType>& mat, unsigned int from, unsigned int to);

            template<typename ElementType>
            static Vector<ElementType> getCol(const Matrix<ElementType>& mat, unsigned int col);

            template<typename ElementType>
            static SparseVector<ElementType> getCol(const SparseMatrix<ElementType>& mat, unsigned int col);

            template<typename ElementType>
            static Vector<ElementType> getCol(const Matrix<ElementType>& mat, unsigned int col, unsigned int from, unsigned int to);

            template<typename ElementType>
            static SparseVector<ElementType> getCol(const SparseMatrix<ElementType>& mat, unsigned int col, unsigned int from, unsigned int to);

            template<typename ElementType>
            static Matrix<ElementType> getCols(const Matrix<ElementType>& mat, unsigned int from, unsigned int to);

            template<typename ElementType>
            static SparseMatrix<ElementType> getCols(const SparseMatrix<ElementType>& mat, unsigned int from, unsigned int to);

            template<typename ElementType>
            static Matrix<ElementType> getSubMatrix(const Matrix<ElementType>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol);

            template<typename ElementType>
            static SparseMatrix<ElementType> getSubMatrix(const SparseMatrix<ElementType>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol);

        };

    }
}

// Include implementation
#include "linalg_traits.hpp"

#endif // GENUSYS_LINALG_LINALG_TRAITS_H_
