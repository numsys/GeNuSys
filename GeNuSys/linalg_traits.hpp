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

#include <algorithm>
#include <complex>

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename SourceType, typename ConvertedType>
        Matrix<ConvertedType> Traits::convertUnsafe(const Matrix<SourceType>& mat)
        {
            Matrix<ConvertedType> result(mat.rows, mat.cols, 00);
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = ElementTraits<SourceType>::template asTypeUnsafe<ConvertedType>(mat.elem[i]);
            }

            return result;
        }

        template<typename SourceType, typename ConvertedType>
        SparseMatrix<ConvertedType> Traits::convertUnsafe(const SparseMatrix<SourceType>& mat)
        {
            SparseMatrix<ConvertedType> result(mat.rows, mat.cols);
            result.row_ptr = mat.row_ptr;
            result.elem = std::vector<typename SparseMatrix<ConvertedType>::Entry>(mat.elem.size());
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                result.elem[i] = typename SparseMatrix<ConvertedType>::Entry(mat.elem[i].col_idx, ElementTraits<SourceType>::template asTypeUnsafe<ConvertedType>(mat.elem[i].value));
            }

            return result;
        }

        template<typename SourceType, typename ConvertedType>
        Vector<ConvertedType> Traits::convertUnsafe(const Vector<SourceType>& vct)
        {
            Vector<ConvertedType> result(vct.length, 00);
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                result.elem[i] = ElementTraits<SourceType>::template asTypeUnsafe<ConvertedType>(vct.elem[i]);
            }

            return result;
        }

        template<typename SourceType, typename ConvertedType>
        SparseVector<ConvertedType> Traits::convertUnsafe(const SparseVector<SourceType>& vct)
        {
            SparseVector<ConvertedType> result(vct.length);
            result.elem = std::vector<typename SparseVector<ConvertedType>::Entry>(vct.elem.size());
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[i] = typename SparseVector<ConvertedType>::Entry(vct.elem[i].idx, ElementTraits<SourceType>::template asTypeUnsafe<ConvertedType>(vct.elem[i].value));
            }

            return result;
        }

        template<typename ElementType>
        void Traits::swapRows(Matrix<ElementType>& mat, const int rowA, const int rowB)
        {
            for (unsigned int i = 0, idxA = rowA * mat.cols, idxB = rowB * mat.cols; i < mat.cols; ++i, ++idxA, ++idxB)
            {
                std::swap(mat.elem[idxA], mat.elem[idxB]);
            }
        }

        template<typename ElementType>
        void Traits::swapCols(Matrix<ElementType>& mat, const int colA, const int colB)
        {
            for (unsigned int i = 0, idxA = colA, idxB = colB; i < mat.rows; ++i, idxA += mat.cols, idxB += mat.cols)
            {
                std::swap(mat.elem[idxA], mat.elem[idxB]);
            }
        }

        template<typename ElementType>
        Vector<ElementType> Traits::getSubVector(const Vector<ElementType>& vct, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= vct.length, std::out_of_range);

            Vector<ElementType> result(to - from, 00);
            for (unsigned int i = 0, idx = from; idx < to; ++i, ++idx)
            {
                result.elem[i] = vct.elem[idx];
            }

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> Traits::getSubVector(const SparseVector<ElementType>& vct, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= vct.length, std::out_of_range);

            SparseVector<ElementType> result(to - from);
            for (unsigned int i = vct.search(from); vct.elem[i].idx < to; ++i)
            {
                result.push(i - from, vct.elem[i].value);
            }

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Traits::getRow(const Matrix<ElementType>& mat, unsigned int row)
        {
            ASSERT_EXCEPTION(row < mat.rows, std::out_of_range);

            Vector<ElementType> result(mat.cols, 00);
            for (unsigned int i = 0, idx = row * mat.cols; i < mat.cols; ++i, ++idx)
            {
                result.elem[i] = mat.elem[idx];
            }

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> Traits::getRow(const SparseMatrix<ElementType>& mat, unsigned int row)
        {
            ASSERT_EXCEPTION(row < mat.rows, std::out_of_range);

            SparseVector<ElementType> result(mat.cols);
            for (unsigned int i = mat.row_ptr[row]; i < mat.row_ptr[row + 1]; ++i)
            {
                result.push(mat.elem[i].col_idx, mat.elem[i].value);
            }

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Traits::getRow(const Matrix<ElementType>& mat, unsigned int row, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(row < mat.rows, std::out_of_range);
            ASSERT_EXCEPTION(from < to && to <= mat.cols, std::out_of_range);

            Vector<ElementType> result(to - from, 00);
            for (unsigned int i = from, idxA = row * mat.cols + from, idxR = 0; i < to; ++i, ++idxA, ++idxR)
            {
                result.elem[idxR] = mat.elem[idxA];
            }

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> Traits::getRow(const SparseMatrix<ElementType>& mat, unsigned int row, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(row < mat.rows, std::out_of_range);
            ASSERT_EXCEPTION(from < to && to <= mat.cols, std::out_of_range);

            SparseVector<ElementType> result(to - from);
            for (unsigned int idx = mat.search(row, from); idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx < to; ++idx)
            {
                result.push(mat.elem[idx].col_idx - from, mat.elem[idx].value);
            }

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Traits::getRows(const Matrix<ElementType>& mat, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= mat.rows, std::out_of_range);

            Matrix<ElementType> result(to - from, mat.cols, 00);
            for (unsigned int i = from, idxA = from * mat.cols, idxR = 0; i < to; ++i)
            {
                for (unsigned int j = 0; j < mat.cols; ++j, ++idxA, ++idxR)
                {
                    result.elem[idxR] = mat.elem[idxA];
                }
            }

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> Traits::getRows(const SparseMatrix<ElementType>& mat, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= mat.rows, std::out_of_range);

            SparseMatrix<ElementType> result(to - from, mat.cols);
            for (unsigned int i = from, rowR = 0; i < to; ++i, ++rowR)
            {
                for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; ++j)
                {
                    result.push(mat.elem[j].idx, mat.elem[j].value);
                }
                result.row_ptr[rowR + 1] = result.size();
            }

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Traits::getCol(const Matrix<ElementType>& mat, unsigned int col)
        {
            ASSERT_EXCEPTION(col < mat.cols, std::out_of_range);

            Vector<ElementType> result(mat.rows, 00);
            for (unsigned int i = 0, idxA = col; i < mat.rows; ++i, idxA += mat.cols)
            {
                result.elem[i] = mat.elem[idxA];
            }
        }

        template<typename ElementType>
        SparseVector<ElementType> Traits::getCol(const SparseMatrix<ElementType>& mat, unsigned int col)
        {
            ASSERT_EXCEPTION(col < mat.cols, std::out_of_range);

            SparseVector<ElementType> result(mat.rows);
            for (unsigned int row = 0; row < mat.rows; ++row)
            {
                unsigned int idx = mat.search(row, col);
                if (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx == col)
                {
                    result.push(row, mat.elem[idx].value);
                }
            }

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Traits::getCol(const Matrix<ElementType>& mat, unsigned int col, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= mat.rows, std::out_of_range);
            ASSERT_EXCEPTION(col < mat.cols, std::out_of_range);

            Vector<ElementType> result(to - from, 00);
            for (unsigned int i = from, idxA = from * mat.cols + col, idxR = 0; i < to; ++i, idxA += mat.cols, ++idxR)
            {
                result.elem[idxR] = mat.elem[idxA];
            }

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> Traits::getCol(const SparseMatrix<ElementType>& mat, unsigned int col, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= mat.rows, std::out_of_range);
            ASSERT_EXCEPTION(col < mat.cols, std::out_of_range);

            SparseVector<ElementType> result(to - from);
            for (unsigned int row = from, rowR = 0; row < to; ++row, ++rowR)
            {
                unsigned int idx = mat.search(row, col);
                if (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx == col)
                {
                    result.push(rowR, mat.elem[idx].value);
                }
                result.row_ptr[rowR + 1] = result.size();
            }

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Traits::getCols(const Matrix<ElementType>& mat, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= mat.cols, std::out_of_range);

            Matrix<ElementType> result(mat.rows, to - from, 00);
            for (unsigned int i = 0, idxA = from, idxR = 0; i < mat.rows; ++i, idxA += mat.cols - to + from)
            {
                for (unsigned int j = from; j < to; ++j, ++idxA, ++idxR)
                {
                    result.elem[idxR] = mat.elem[idxA];
                }
            }

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> Traits::getCols(const SparseMatrix<ElementType>& mat, unsigned int from, unsigned int to)
        {
            ASSERT_EXCEPTION(from < to && to <= mat.cols, std::out_of_range);

            SparseMatrix<ElementType> result(mat.rows, to - from);
            for (unsigned int row = 0; row < mat.rows; ++row)
            {
                unsigned int idx = mat.search(row, from);
                while (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx < to)
                {
                    result.push(mat.elem[idx].col_idx - from, mat.elem[idx].value);
                }
                result.row_ptr[row + 1] = result.size();
            }

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Traits::getSubMatrix(const Matrix<ElementType>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol)
        {
            ASSERT_EXCEPTION(fromRow < toRow && toRow <= mat.rows, std::out_of_range);
            ASSERT_EXCEPTION(fromCol < toCol && toCol <= mat.cols, std::out_of_range);

            Matrix<ElementType> result(toRow - fromRow, toCol - fromCol, 00);
            for (unsigned int i = fromRow, idxA = fromRow * mat.cols + fromCol, idxR = 0; i < toRow; ++i, idxA += mat.cols - toCol + fromCol)
            {
                for (unsigned int j = fromCol; j < toCol; ++j, ++idxA, ++idxR)
                {
                    result.elem[idxR] = mat.elem[idxA];
                }
            }

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> Traits::getSubMatrix(const SparseMatrix<ElementType>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol)
        {
            ASSERT_EXCEPTION(fromRow < toRow && toRow <= mat.rows, std::out_of_range);
            ASSERT_EXCEPTION(fromCol < toCol && toCol <= mat.cols, std::out_of_range);

            SparseMatrix<ElementType> result(toRow - fromRow, toCol - fromCol);
            for (unsigned int row = fromRow, rowR = 0; row < toRow; ++row, ++rowR)
            {
                unsigned int idx = mat.search(row, fromCol);
                while (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx < toCol)
                {
                    result.push(mat.elem[idx].col_idx - fromCol, mat.elem[idx].value);
                }
                result.row_ptr[rowR + 1] = result.size();
            }

            return result;
        }

    }
}
