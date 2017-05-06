#include <stdexcept>
#include <algorithm>

#include "linalg_operations.h"
#include "utils.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::diag(const Vector<ElementType>& vct)
        {
            Matrix<ElementType> result(vct.length, vct.length);
            for (unsigned int i = 0, idx = 0; i < vct.length; ++i, idx += vct.length + 1)
            {
                result.elem[idx] = vct.elem[i];
            }

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::diag(const SparseVector<ElementType>& vct)
        {
            Matrix<ElementType> result(vct.length, vct.length);
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                result.elem[vct.elem[i].idx * (vct.length + 1)] = vct.elem[i].value;
            }

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::identity(unsigned int rows, unsigned int cols)
        {
            Matrix<ElementType> result(rows, cols);
            unsigned int min = cols < rows ? cols : rows;
            for (unsigned int i = 0, idx = 0; i < min; ++i, idx += cols + 1)
            {
                result.elem[idx] = ElementTraits<ElementType>::one();
            }

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType>::Matrix(): rows(0), cols(0), elem(rows * cols)
        {
        }

        template<typename ElementType>
        Matrix<ElementType>::Matrix(unsigned int rows, unsigned int cols): rows(rows), cols(cols), elem(rows * cols, ElementTraits<ElementType>::zero())
        {
        }

        template<typename ElementType>
        Matrix<ElementType>::Matrix(unsigned int rows, unsigned int cols, int): rows(rows), cols(cols), elem(rows * cols)
        {
        }

        template<typename ElementType>
        Matrix<ElementType>::Matrix(const Matrix<ElementType>& mat): rows(mat.rows), cols(mat.cols), elem(mat.elem)
        {
        }

        template<typename ElementType>
        template<typename SourceType>
        Matrix<ElementType>::Matrix(const Matrix<SourceType>& mat): rows(mat.rows), cols(mat.cols), elem(mat.rows * mat.cols)
        {
            for (unsigned int i = 0; i < size(); ++i)
            {
                elem[i] = ElementTraits<SourceType>::template asType<ElementType>(mat.elem[i]);
            }
        }

        template<typename ElementType>
        template<typename SourceType>
        Matrix<ElementType>::Matrix(const SparseMatrix<SourceType>& mat): rows(mat.rows), cols(mat.cols), elem(mat.rows * mat.cols, ElementTraits<SourceType>::zero())
        {
            for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i, idxA += cols)
            {
                for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; ++j)
                {
                    elem[idxA + mat.elem[j].col_idx] = ElementTraits<SourceType>::template asType<ElementType>(mat.elem[j].value);
                }
            }
        }

        template<typename ElementType>
        Matrix<ElementType>::~Matrix()
        {
        }

        template<typename ElementType>
        Matrix<ElementType>& Matrix<ElementType>::operator =(const Matrix<ElementType>& mat)
        {
            if (this != &mat)
            {
                rows = mat.rows;
                cols = mat.cols;
                elem = mat.elem;
            }

            return *this;
        }

        template<typename ElementType>
        template<typename SourceType>
        Matrix<ElementType>& Matrix<ElementType>::operator =(const Matrix<SourceType>& mat)
        {
            rows = mat.rows;
            cols = mat.cols;
            elem = std::vector<ElementType>(mat.rows * mat.cols);
            for (unsigned int i = 0; i < size(); ++i)
            {
                elem[i] = ElementTraits<SourceType>::template asType<ElementType>(mat.elem[i]);
            }

            return *this;
        }

        template<typename ElementType>
        template<typename SourceType>
        Matrix<ElementType>& Matrix<ElementType>::operator =(const SparseMatrix<SourceType>& mat)
        {
            rows = mat.rows;
            cols = mat.cols;
            elem = std::vector<ElementType>(mat.rows * mat.cols, ElementTraits<SourceType>::zero());
            for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i, idxA += cols)
            {
                for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; ++j)
                {
                    elem[idxA + mat.elem[j].col_idx] = ElementTraits<SourceType>::template asType<ElementType>(mat.elem[j].value);
                }
            }

            return *this;
        }

        template<typename ElementType>
        unsigned int Matrix<ElementType>::size() const
        {
            return elem.size();
        }

        template<typename ElementType>
        unsigned int Matrix<ElementType>::getRows() const
        {
            return rows;
        }

        template<typename ElementType>
        unsigned int Matrix<ElementType>::getCols() const
        {
            return cols;
        }

        template<typename ElementType>
        void Matrix<ElementType>::set(unsigned int i, unsigned int j, const ElementType& value)
        {
            ASSERT_EXCEPTION(i < rows && j < cols, std::out_of_range);

            elem[i * cols + j] = value;
        }

        template<typename ElementType>
        const ElementType& Matrix<ElementType>::operator()(unsigned int i, unsigned int j) const
        {
            ASSERT_EXCEPTION(i < rows && j < cols, std::out_of_range);

            return elem[i * cols + j];
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::transpose() const
        {
            Matrix<ElementType> result(cols, rows, 00);
            Operations::mat_transpose(*this, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::conjugateTranspose() const
        {
            Matrix<ElementType> result(cols, rows, 00);
            Operations::mat_conjugate_transpose(*this, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::operator *(const ElementType& value) const
        {
            Matrix<ElementType> result(rows, cols, 00);
            Operations::mat_mul(*this, value, result);

            return result;
        }

        template<typename ElementType>
        Matrix<typename ElementTraits<ElementType>::RationalType> Matrix<ElementType>::operator /(const ElementType& value) const
        {
            Matrix<typename ElementTraits<ElementType>::RationalType> result(rows, cols, 00);
            Operations::mat_div(*this, value, result);

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Matrix<ElementType>::operator *(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            Vector<ElementType> result(rows, 00);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> Matrix<ElementType>::operator *=(const Vector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            SparseVector<ElementType> result(rows);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        Vector<ElementType> Matrix<ElementType>::operator *(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            Vector<ElementType> result(rows, 00);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        SparseVector<ElementType> Matrix<ElementType>::operator *=(const SparseVector<ElementType>& vct) const
        {
            ASSERT_EXCEPTION(cols == vct.length, std::length_error);

            SparseVector<ElementType> result(rows);
            Operations::mat_mul(*this, vct, result);

            return result;
        }

        template<typename ElementType>
        bool Matrix<ElementType>::operator ==(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return Operations::mat_eq(*this, mat);
        }

        template<typename ElementType>
        bool Matrix<ElementType>::operator !=(const Matrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return !Operations::mat_eq(*this, mat);
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::operator +(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            Matrix<ElementType> result(rows, cols, 00);
            Operations::mat_add(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::operator -(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            Matrix<ElementType> result(rows, cols, 00);
            Operations::mat_sub(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::operator *(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            Matrix<ElementType> result(rows, mat.cols, 00);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> Matrix<ElementType>::operator *=(const Matrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            SparseMatrix<ElementType> result(rows, mat.cols);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::operator +(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            Matrix<ElementType> result(rows, cols);
            Operations::mat_add(*this, mat, result);


            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::operator -(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            Matrix<ElementType> result(rows, cols, 00);
            Operations::mat_sub(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        Matrix<ElementType> Matrix<ElementType>::operator *(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            Matrix<ElementType> result(rows, mat.cols, 00);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        SparseMatrix<ElementType> Matrix<ElementType>::operator *=(const SparseMatrix<ElementType>& mat) const
        {
            ASSERT_EXCEPTION(cols == mat.rows, std::length_error);

            SparseMatrix<ElementType> result(rows, mat.cols);
            Operations::mat_mul(*this, mat, result);

            return result;
        }

        template<typename ElementType>
        bool Matrix<ElementType>::operator ==(const SparseMatrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return Operations::mat_eq(*this, mat);
        }

        template<typename ElementType>
        bool Matrix<ElementType>::operator !=(const SparseMatrix<ElementType>& mat)
        {
            ASSERT_EXCEPTION(rows == mat.rows && cols == mat.cols, std::length_error);

            return !Operations::mat_eq(*this, mat);
        }

    }
}
