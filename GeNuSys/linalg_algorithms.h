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

#ifndef GENUSYS_LINALG_LINALG_ALGORITHMS_H_
#define GENUSYS_LINALG_LINALG_ALGORITHMS_H_

#include <vector>

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<typename ElementType>
        struct LU
        {

            const Matrix<typename ElementTraits<ElementType>::RationalType> L;

            const Matrix<typename ElementTraits<ElementType>::RationalType> U;

            LU(const Matrix<typename ElementTraits<ElementType>::RationalType>& L,
               const Matrix<typename ElementTraits<ElementType>::RationalType>& U);

        };

        template<typename ElementType>
        struct QR
        {

            const Matrix<typename ElementTraits<ElementType>::RealType> Q;

            const Matrix<typename ElementTraits<ElementType>::RealType> R;

            QR(const Matrix<typename ElementTraits<ElementType>::RealType>& Q,
               const Matrix<typename ElementTraits<ElementType>::RealType>& R);

        };

        template<typename ElementType>
        struct SmithNormalForm
        {

            const Matrix<ElementType> S;

            const Matrix<ElementType> U;

            const Matrix<ElementType> V;

            SmithNormalForm(const Matrix<ElementType>& S,
                            const Matrix<ElementType>& U,
                            const Matrix<ElementType>& V);

        };

        template<typename ElementType>
        struct HessenbergForm
        {

            const Matrix<typename ElementTraits<ElementType>::RealType> Q;

            const Matrix<typename ElementTraits<ElementType>::RealType> A;

            const Matrix<typename ElementTraits<ElementType>::RealType> QT;

            HessenbergForm(const Matrix<typename ElementTraits<ElementType>::RealType>& Q,
                           const Matrix<typename ElementTraits<ElementType>::RealType>& A,
                           const Matrix<typename ElementTraits<ElementType>::RealType>& QT);

        };

        template<typename ElementType>
        struct SchurForm
        {

            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType> Q;

            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType> U;

            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType> QT;

            SchurForm(const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& Q,
                      const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& U,
                      const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& QT);

        };

        template<typename ElementType>
        struct JordanForm
        {

            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType> P;

            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType> J;

            const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType> invP;

            JordanForm(const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& P,
                       const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& J,
                       const Matrix<typename ElementTraits<typename ElementTraits<ElementType>::ComplexType>::RealType>& invP);

        };

        struct Algorithms
        {

            template<typename ElementType>
            static typename ElementTraits<ElementType>::RationalType det(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static Matrix<typename ElementTraits<ElementType>::RationalType> invert(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static Matrix<ElementType> getAdjoint(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static LU<ElementType> decomposeLU(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static Matrix<typename ElementTraits<ElementType>::RealType> getHouseholderMatrix(const Vector<ElementType>& vct, const unsigned int N);

            template<typename ElementType>
            static void decomposeQR(const Matrix<ElementType>& mat, QR<ElementType>& qr);

            template<typename ElementType>
            static QR<ElementType> decomposeQR(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static SmithNormalForm<ElementType> getSmithNormalForm(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static HessenbergForm<ElementType> getHessenbergForm(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static SchurForm<ElementType> getSchurForm(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static unsigned int getRank(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static std::vector<Vector<typename ElementTraits<ElementType>::RationalType>> solveHomogeneous(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static std::vector<Vector<typename ElementTraits<ElementType>::RationalType>> getEigenVectors(const Matrix<ElementType>& mat, const ElementType& eigenValue);

            template<typename ElementType>
            static JordanForm<ElementType> getJordanForm(const Matrix<ElementType>& mat);

        };

    }
}

// Include implementation
#include "linalg_algorithms.hpp"

#endif // GENUSYS_LINALG_LINALG_ALGORITHMS_H_
