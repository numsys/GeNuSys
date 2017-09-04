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

#ifndef GENUSYS_LINALG_OPERATOR_NORM_H_
#define GENUSYS_LINALG_OPERATOR_NORM_H_

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

#include "linalg_algorithms.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename BaseType>
        class OperatorNorm
        {

            private:

                typedef typename ElementTraits<typename ElementTraits<BaseType>::ComplexType>::RealType ComplexRealType;

                Matrix<ComplexRealType> S;
                Matrix<ComplexRealType> invS;

            public:

                OperatorNorm(const Matrix<BaseType>& mat);

                OperatorNorm(const JordanForm<BaseType>& jordanForm);

                template<typename ElementType>
                struct NormType
                {

                    typedef typename PNorm<00>::template NormType<ComplexRealType>::Type Type;

                };

                template<typename ElementType>
                typename NormType<ElementType>::Type norm(const Vector<ElementType>& vct) const;

                template<typename ElementType>
                typename NormType<ElementType>::Type norm(const SparseVector<ElementType>& vct) const;

                template<typename ElementType>
                typename NormType<ElementType>::Type norm(const Matrix<ElementType>& mat) const;

                template<typename ElementType>
                typename NormType<ElementType>::Type norm(const SparseMatrix<ElementType>& mat) const;

                template<typename ElementType>
                typename NormType<ElementType>::Type operator()(const Vector<ElementType>& vct) const;

                template<typename ElementType>
                typename NormType<ElementType>::Type operator()(const SparseVector<ElementType>& vct) const;

                template<typename ElementType>
                typename NormType<ElementType>::Type operator()(const Matrix<ElementType>& mat) const;

                template<typename ElementType>
                typename NormType<ElementType>::Type operator()(const SparseMatrix<ElementType>& mat) const;

        };

    }
}

// Include implementation
#include "operator_norm.hpp"

#endif // GENUSYS_LINALG_OPERATOR_NORM_H_
