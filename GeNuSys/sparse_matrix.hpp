#include <cassert>
#include <algorithm>

#include "linalg_operations.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        SparseMatrix<ElementType>::Entry::Entry()
        {
        }

        template<typename ElementType>
        SparseMatrix<ElementType>::Entry::Entry(unsigned int col_idx, const ElementType& value): col_idx(col_idx), value(value)
        {
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::diag(const Vector<ElementType>& vct)
        {
            SparseMatrix<ElementType> result(vct.length, vct.length);
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                result.elem[i] = Entry(i, vct.elem[i]);
                result.row_ptr[i + 1] = i + 1;
            }

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::diag(const SparseVector<ElementType>& vct)
        {
            SparseMatrix<ElementType> result(vct.length, vct.length, vct.size());
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[i] = Entry(vct.elem[i].idx, vct.elem[i].elem);
                result.row_ptr[i + 1] = i + 1;
            }

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::identity(unsigned int rows, unsigned int cols)
        {
            unsigned int min = cols < rows ? cols : rows;

            SparseMatrix<ElementType> result(rows, cols, min);
            for (unsigned int i = 0; i < min; ++i)
            {
                result.elem[i] = Entry(i, ElementTraits<ElementType>::one());
                result.row_ptr[i + 1] = i + 1;
            }
            for (unsigned int i = min; i < rows; ++i)
            {
                result.row_ptr[i + 1] = i;
            }

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType>::SparseMatrix(): rows(0), cols(0), row_ptr(0), elem()
        {
        }

        template<typename ElementType>
        SparseMatrix<ElementType>::SparseMatrix(unsigned int rows, unsigned int cols, unsigned int capacity): rows(rows), cols(cols), row_ptr(rows + 1, 0), elem(capacity)
        {
        }

        template<typename ElementType>
        SparseMatrix<ElementType>::SparseMatrix(unsigned int rows, unsigned int cols): rows(rows), cols(cols), row_ptr(rows + 1, 0), elem()
        {
        }

        template<typename ElementType>
        SparseMatrix<ElementType>::SparseMatrix(const SparseMatrix<ElementType>& mat): rows(mat.rows), cols(mat.cols), row_ptr(mat.row_ptr), elem(mat.elem)
        {
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseMatrix<ElementType>::SparseMatrix(const SparseMatrix<SourceType>& mat): rows(mat.rows), cols(mat.cols), row_ptr(mat.row_ptr), elem(mat.elem.size())
        {
            for (unsigned int i = 0; i < size(); ++i)
            {
                elem[i] = Entry(mat.elem[i].col_idx, ElementTraits<SourceType>::template asType<ElementType>(mat.elem[i].value));
            }
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseMatrix<ElementType>::SparseMatrix(const Matrix<SourceType>& mat): rows(mat.rows), cols(mat.cols), row_ptr(rows + 1, 0)
        {
            unsigned int nnz = 0;
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                if (mat.elem[i] != ElementTraits<SourceType>::zero())
                {
                    ++nnz;
                }
            }

            elem = std::vector<Entry>(nnz);
            unsigned int size = 0;
            for (unsigned int i = 0, idxB = 0; i < rows; ++i)
            {
                for (unsigned int j = 0; j < cols; ++j, ++idxB)
                {
                    if (mat.elem[idxB] != ElementTraits<SourceType>::zero())
                    {
                        elem[size++] = Entry(j, ElementTraits<SourceType>::template asType<ElementType>(mat.elem[idxB]));
                    }
                }
                row_ptr[i + 1] = size;
            }
        }

        template<typename ElementType>
        SparseMatrix<ElementType>::~SparseMatrix()
        {
        }

        template<typename ElementType>
        SparseMatrix<ElementType>& SparseMatrix<ElementType>::operator =(const SparseMatrix<ElementType>& mat)
        {
            if (this != &mat)
            {
                rows = mat.rows;
                cols = mat.cols;
                row_ptr = mat.row_ptr;
                elem = mat.elem;
            }

            return *this;
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseMatrix<ElementType>& SparseMatrix<ElementType>::operator =(const SparseMatrix<SourceType>& mat)
        {
            rows = mat.rows;
            cols = mat.cols;
            row_ptr = mat.row_ptr;
            elem = std::vector<Entry>(mat.elem.size());
            for (unsigned int i = 0; i < size(); ++i)
            {
                elem[i] = Entry(mat.elem[i].col_idx, ElementTraits<SourceType>::template asType<ElementType>(mat.elem[i].value));
            }

            return *this;
        }

        template<typename ElementType>
        template<typename SourceType>
        SparseMatrix<ElementType>& SparseMatrix<ElementType>::operator =(const Matrix<SourceType>& mat)
        {
            rows = mat.rows;
            cols = mat.cols;
            row_ptr = std::vector<unsigned int>(rows + 1, 0);

            unsigned int nnz = 0;
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                if (mat.elem[i] != ElementTraits<SourceType>::zero())
                {
                    ++nnz;
                }
            }

            elem = std::vector<Entry>(nnz);
            unsigned int size = 0;
            for (unsigned int i = 0, idxB = 0; i < rows; ++i)
            {
                for (unsigned int j = 0; j < cols; ++j, ++idxB)
                {
                    if (mat.elem[idxB] != ElementTraits<SourceType>::zero())
                    {
                        elem[size++] = Entry(j, ElementTraits<SourceType>::template asType<ElementType>(mat.elem[idxB]));
                    }
                }
                row_ptr[i + 1] = size;
            }

            return *this;
        }

        template<typename ElementType>
        unsigned int SparseMatrix<ElementType>::size() const
        {
            return elem.size();
        }

        template<typename ElementType>
        unsigned int SparseMatrix<ElementType>::search(unsigned int row_idx, unsigned int col_idx) const
        {
            unsigned int min = row_ptr[row_idx];
            unsigned int max = row_ptr[row_idx + 1];
            unsigned int mid;
            while (min < max)
            {
                mid = (min + max) / 2;
                if (elem[mid].col_idx < col_idx)
                {
                    min = mid + 1;
                }
                else
                {
                    max = mid;
                }
            }
            return min;
        }

        template<typename ElementType>
        void SparseMatrix<ElementType>::push(unsigned int col_idx, const ElementType& value)
        {
            elem.push_back(Entry(col_idx, value));
        }

        template<typename ElementType>
        unsigned int SparseMatrix<ElementType>::getRows() const
        {
            return rows;
        }

        template<typename ElementType>
        unsigned int SparseMatrix<ElementType>::getCols() const
        {
            return cols;
        }

        template<typename ElementType>
        void SparseMatrix<ElementType>::set(unsigned int row_idx, unsigned int col_idx, const ElementType& value)
        {
            ASSERT_EXCEPTION(row_idx < rows && col_idx < cols, std::out_of_range);

            unsigned int entry_idx = search(row_idx, col_idx);
            if (value != ElementTraits<ElementType>::zero())
            {
                if (entry_idx < row_ptr[row_idx + 1] && elem[entry_idx].col_idx == col_idx)
                {
                    elem[entry_idx].value = value;
                }
                else
                {
                    elem.insert(elem.begin() + entry_idx, Entry(col_idx, value));
                    for (unsigned int i = row_idx + 1; i <= rows; ++i)
                    {
                        ++row_ptr[i];
                    }
                }
            }
            else
            {
                if (entry_idx < row_ptr[row_idx + 1] && elem[entry_idx].col_idx == col_idx)
                {
                    elem.erase(elem.begin() + entry_idx);
                    for (unsigned int i = row_idx + 1; i <= rows; ++i)
                    {
                        ++row_ptr[i];
                    }
                }
            }
        }

        template<typename ElementType>
        const ElementType& SparseMatrix<ElementType>::operator()(unsigned int row_idx, unsigned int col_idx) const
        {
            ASSERT_EXCEPTION(row_idx < rows && col_idx < cols, std::out_of_range);

            unsigned int entry_idx = search(row_idx, col_idx);
            if (entry_idx < row_ptr[row_idx + 1] && elem[entry_idx].col_idx == col_idx)
            {
                return elem[entry_idx].value;
            }
            else
            {
                return ElementTraits<ElementType>::zero();
            }
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::transpose() const
        {
            SparseMatrix<ElementType> result(cols, rows, size());
            Operations::mat_transpose(*this, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::conjugateTranspose() const
        {
            SparseMatrix<ElementType> result(cols, rows, size());
            Operations::mat_conjugate_transpose(*this, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::operator *(const ElementType& value) const
        {
            SparseMatrix<ElementType> result(rows, cols, size());
            Operations::mat_mul(*this, value, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<typename ElementTraits<ElementType>::RationalType> SparseMatrix<ElementType>::operator /(const ElementType& value) const
        {
            typedef typename ElementTraits<ElementType>::RationalType RationalType;

            SparseMatrix<RationalType> result(rows, cols, size());
            Operations::mat_div(*this, value, result);

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> SparseMatrix<ElementType>::operator *(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            Vector<ElementType> result(rows, 00);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> SparseMatrix<ElementType>::operator *=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            SparseVector<ElementType> result(rows);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> SparseMatrix<ElementType>::operator *(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            Vector<ElementType> result(rows, 00);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> SparseMatrix<ElementType>::operator *=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            SparseVector<ElementType> result(rows);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> SparseMatrix<ElementType>::operator +(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            Matrix<ElementType> result(rows, cols, 00);
            Operations::mat_add(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> SparseMatrix<ElementType>::operator -(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            Matrix<ElementType> result(rows, cols, 00);
            Operations::mat_sub(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> SparseMatrix<ElementType>::operator *(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            Matrix<ElementType> result(rows, mat.cols, 00);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::operator *=(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            SparseMatrix<ElementType> result(rows, mat.cols);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        bool SparseMatrix<ElementType>::operator ==(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return Operations::mat_eq(*this, mat);;
        }

        template<typename ElementType>
        bool SparseMatrix<ElementType>::operator !=(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return !Operations::mat_eq(*this, mat);
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::operator +(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            SparseMatrix<ElementType> result(rows, cols);
            Operations::mat_add(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::operator -(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            SparseMatrix<ElementType> result(rows, cols);
            Operations::mat_sub(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> SparseMatrix<ElementType>::operator *(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            Matrix<ElementType> result(rows, mat.cols, 00);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> SparseMatrix<ElementType>::operator *=(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            SparseMatrix<ElementType> result(rows, mat.cols);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        bool SparseMatrix<ElementType>::operator ==(const SparseMatrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return Operations::mat_eq(*this, mat);
        }

        template<typename ElementType>
        bool SparseMatrix<ElementType>::operator !=(const SparseMatrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return !Operations::mat_eq(*this, mat);
        }

    }
}
