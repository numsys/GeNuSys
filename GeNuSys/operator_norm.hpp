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

#include "linalg_algorithms.h"
#include "p_norm.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename BaseType>
        OperatorNorm<BaseType>::OperatorNorm(const Matrix<BaseType>& mat)
        {
            typedef typename ElementTraits<ComplexRealType>::AbsType AbsType;

            JordanForm<BaseType> jordanForm = Algorithms::getJordanForm(mat);

            const unsigned int N = jordanForm.J.rows;

            Matrix<ComplexRealType> D = Matrix<ComplexRealType>::identity(N, N);

            bool first = true;
            ComplexRealType act;
            AbsType mu, muPow;
            for (int i = N - 1; i >= 0; ++i)
            {
                if (first || jordanForm.J(i, i) != act)
                {
                    mu = ElementTraits<AbsType>::one() - ElementTraits<ComplexRealType>::abs(jordanForm.J(i, i));
                    muPow = mu;
                    first = false;
                }
                else
                {
                    muPow *= mu;
                }
                D.set(i, i, muPow);
            }

            S = jordanForm.P * D;
            invS = Algorithms::invert(S);
        }

        template<typename BaseType>
        OperatorNorm<BaseType>::OperatorNorm(const JordanForm<BaseType>& jordanForm)
        {
            typedef typename ElementTraits<ComplexRealType>::AbsType AbsType;

            const unsigned int N = jordanForm.J.rows;

            Matrix<ComplexRealType> D = Matrix<ComplexRealType>::identity(N, N);

            bool first = true;
            ComplexRealType act;
            AbsType mu, muPow;
            for (int i = N - 1; i >= 0; --i)
            {
                if (first || jordanForm.J(i, i) != act)
                {
                    mu = ElementTraits<AbsType>::one() - ElementTraits<ComplexRealType>::abs(jordanForm.J(i, i));
                    if (mu <= ElementTraits<AbsType>::zero())
                    {
                        mu = ElementTraits<AbsType>::one();
                    }
                    muPow = mu;
                    first = false;
                }
                else
                {
                    muPow *= mu;
                }
                D.set(i, i, muPow);
            }

            S = jordanForm.P * D;
            invS = Algorithms::invert(S);
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::norm(const Vector<ElementType>& vct) const
        {
            return PNorm<00>::norm(S * Vector<ComplexRealType>(vct));
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::norm(const SparseVector<ElementType>& vct) const
        {
            return PNorm<00>::norm(S * SparseVector<ComplexRealType>(vct));
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::norm(const Matrix<ElementType>& mat) const
        {
            return PNorm<00>::norm(S * Matrix<ComplexRealType>(mat) * invS);
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::norm(const SparseMatrix<ElementType>& mat) const
        {
            return PNorm<00>::norm(S * SparseMatrix<ComplexRealType>(mat) * invS);
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::operator()(const Vector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::operator()(const SparseVector<ElementType>& vct) const
        {
            return norm(vct);
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::operator()(const Matrix<ElementType>& mat) const
        {
            return norm(mat);
        }

        template<typename BaseType>
        template<typename ElementType>
        typename OperatorNorm<BaseType>::template NormType<ElementType>::Type OperatorNorm<BaseType>::operator()(const SparseMatrix<ElementType>& mat) const
        {
            return norm(mat);
        }

    }
}
