#include <vector>

namespace GeNuSys
{
    namespace LinAlg
    {

        template<unsigned int P>
        template<typename ElementType>
        typename PNorm<P>::template NormType<ElementType>::Type PNorm<P>::norm(const Vector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType sum = ElementTraits<AbsType>::zero();
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                sum += ElementTraits<AbsType>::pow(ElementTraits<ElementType>::abs(vct.elem[i]), P);
            }

            return ElementTraits<AbsType>::root(sum, P);
        }

        template<unsigned int P>
        template<typename ElementType>
        typename PNorm<P>::template NormType<ElementType>::Type PNorm<P>::norm(const SparseVector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType sum = ElementTraits<AbsType>::zero();
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                sum += ElementTraits<AbsType>::pow(ElementTraits<ElementType>::abs(vct.elem[i].value), P);
            }

            return ElementTraits<AbsType>::root(sum, P);
        }

        template<unsigned int P>
        template<typename ElementType>
        typename PNorm<P>::template NormType<ElementType>::Type PNorm<P>::operator()(const Vector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<unsigned int P>
        template<typename ElementType>
        typename PNorm<P>::template NormType<ElementType>::Type PNorm<P>::operator()(const SparseVector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<>
        struct PNorm<1>
        {

            template<typename ElementType>
            struct NormType
            {

                typedef typename ElementTraits<ElementType>::AbsType Type;

            };

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Vector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseVector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseMatrix<ElementType>& mat);

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Vector<ElementType>& vct) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseVector<ElementType>& vct) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Matrix<ElementType>& mat) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseMatrix<ElementType>& mat) const;

        };

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::norm(const Vector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType sum = ElementTraits<AbsType>::zero();
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                sum += ElementTraits<ElementType>::abs(vct.elem[i]);
            }

            return sum;
        }

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::norm(const SparseVector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType sum = ElementTraits<AbsType>::zero();
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                sum += ElementTraits<ElementType>::abs(vct.elem[i].value);
            }

            return sum;
        }

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::norm(const Matrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType max = ElementTraits<AbsType>::zero();
            for (unsigned int row = 0, idx = 0; row < mat.rows; ++row, idx += mat.cols)
            {
                max += ElementTraits<ElementType>::abs(mat.elem[idx]);
            }
            for (unsigned int col = 1; col < mat.cols; ++col)
            {
                AbsType sum = ElementTraits<AbsType>::zero();
                for (unsigned int row = 0, idx = col; row < mat.rows; ++row, idx += mat.cols)
                {
                    sum += ElementTraits<ElementType>::abs(mat.elem[idx]);
                }
                if (max < sum)
                {
                    max = sum;
                }
            }

            return max;
        }

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::norm(const SparseMatrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            std::vector<AbsType> cache(mat.cols, ElementTraits<AbsType>::zero());
            for (unsigned int i = 0; i < mat.size(); ++i)
            {
                cache[mat.elem[i].col_idx] += ElementTraits<ElementType>::abs(mat.elem[i].value);
            }

            AbsType max = cache[0];
            for (unsigned int i = 1; i < mat.cols; ++i)
            {
                if (max < cache[i])
                {
                    max = cache[i];
                }
            }

            return max;
        }

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::operator()(const Vector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::operator()(const SparseVector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::operator()(const Matrix<ElementType>& mat) const
        {
            return norm(mat);
        }

        template<typename ElementType>
        typename PNorm<1>::template NormType<ElementType>::Type PNorm<1>::operator()(const SparseMatrix<ElementType>& mat) const
        {
            return norm(mat);
        }

        template<>
        struct PNorm<2>
        {

            template<typename ElementType>
            struct NormType
            {

                typedef typename ElementTraits<ElementType>::AbsSqrType SqrType;

                typedef typename ElementTraits<SqrType>::RealType Type;

            };

            template<typename ElementType>
            static typename NormType<ElementType>::SqrType normSqr(const Vector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::SqrType normSqr(const SparseVector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Vector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseVector<ElementType>& vct);

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Vector<ElementType>& vct) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseVector<ElementType>& vct) const;

        };

        template<typename ElementType>
        typename PNorm<2>::template NormType<ElementType>::SqrType PNorm<2>::normSqr(const Vector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsSqrType AbsSqrType;

            AbsSqrType sum = ElementTraits<AbsSqrType>::zero();
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                sum += ElementTraits<ElementType>::absSqr(vct.elem[i]);
            }

            return sum;
        }

        template<typename ElementType>
        typename PNorm<2>::template NormType<ElementType>::SqrType PNorm<2>::normSqr(const SparseVector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsSqrType AbsSqrType;

            AbsSqrType sum = ElementTraits<AbsSqrType>::zero();
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                sum += ElementTraits<ElementType>::absSqr(vct.elem[i].value);
            }

            return sum;
        }

        template<typename ElementType>
        typename PNorm<2>::template NormType<ElementType>::Type PNorm<2>::norm(const Vector<ElementType>& vct)
        {
            return ElementTraits<typename NormType<ElementType>::SqrType>::sqrt(normSqr(vct));
        }

        template<typename ElementType>
        typename PNorm<2>::template NormType<ElementType>::Type PNorm<2>::norm(const SparseVector<ElementType>& vct)
        {
            return ElementTraits<typename NormType<ElementType>::SqrType>::sqrt(normSqr(vct));
        }

        template<typename ElementType>
        typename PNorm<2>::template NormType<ElementType>::Type PNorm<2>::operator()(const Vector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<typename ElementType>
        typename PNorm<2>::template NormType<ElementType>::Type PNorm<2>::operator()(const SparseVector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<>
        struct PNorm<00>
        {

            template<typename ElementType>
            struct NormType
            {

                typedef typename ElementTraits<ElementType>::AbsType Type;

            };

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Vector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseVector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseMatrix<ElementType>& mat);

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Vector<ElementType>& vct) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseVector<ElementType>& vct) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Matrix<ElementType>& mat) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseMatrix<ElementType>& mat) const;

        };

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::norm(const Vector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType max = ElementTraits<AbsType>::zero();
            for (unsigned int i = 0; i < vct.length; ++i)
            {
                AbsType abs = ElementTraits<ElementType>::abs(vct.elem[i]);
                if (abs > max)
                {
                    max = abs;
                }
            }

            return max;
        }

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::norm(const SparseVector<ElementType>& vct)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType max = ElementTraits<AbsType>::zero();
            for (unsigned int i = 0; i < vct.size(); ++i)
            {
                AbsType abs = ElementTraits<ElementType>::abs(vct.elem[i].value);
                if (abs > max)
                {
                    max = abs;
                }
            }

            return max;
        }

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::norm(const Matrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            unsigned int idx = 0;
            AbsType max = ElementTraits<AbsType>::zero();
            for (unsigned int col = 0; col < mat.cols; ++col)
            {
                max += ElementTraits<ElementType>::abs(mat.elem[idx++]);
            }
            for (unsigned int row = 1; row < mat.rows; ++row)
            {
                AbsType sum = ElementTraits<AbsType>::zero();
                for (unsigned int col = 0; col < mat.cols; ++col)
                {
                    sum += ElementTraits<ElementType>::abs(mat.elem[idx++]);
                }
                if (max < sum)
                {
                    max = sum;
                }
            }

            return max;
        }

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::norm(const SparseMatrix<ElementType>& mat)
        {
            typedef typename ElementTraits<ElementType>::AbsType AbsType;

            AbsType max = ElementTraits<AbsType>::zero();
            for (unsigned int idx = mat.row_ptr[0]; idx < mat.row_ptr[1]; ++idx)
            {
                max += ElementTraits<ElementType>::abs(mat.elem[idx].value);
            }
            for (unsigned int row = 1; row < mat.rows; ++row)
            {
                AbsType sum = ElementTraits<AbsType>::zero();
                for (unsigned int idx = mat.row_ptr[row]; idx < mat.row_ptr[row + 1]; ++idx)
                {
                    sum += ElementTraits<ElementType>::abs(mat.elem[idx].value);
                }
                if (max < sum)
                {
                    max = sum;
                }
            }

            return max;
        }

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::operator()(const Vector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::operator()(const SparseVector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::operator()(const Matrix<ElementType>& mat) const
        {
            return norm(mat);
        }

        template<typename ElementType>
        typename PNorm<00>::template NormType<ElementType>::Type PNorm<00>::operator()(const SparseMatrix<ElementType>& mat) const
        {
            return norm(mat);
        }

    }
}
