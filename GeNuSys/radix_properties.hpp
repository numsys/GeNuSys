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

#include "element_traits.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template<typename ElementType>
        RadixProperties<ElementType>::RadixProperties(const GeNuSys::LinAlg::Matrix<ElementType>& M):
            M(M),
            invM(GeNuSys::LinAlg::Algorithms::invert(M)),
            adjM(GeNuSys::LinAlg::Algorithms::getAdjoint(M)),
            detM(ElementTraits<typename ElementTraits<ElementType>::RationalType>::template asTypeUnsafe<ElementType>(GeNuSys::LinAlg::Algorithms::det(M))),
            absDetM(ElementTraits<ElementType>::abs(detM)),
            smithNormalForm(GeNuSys::LinAlg::Algorithms::getSmithNormalForm(M)),
            operatorNorm(GeNuSys::LinAlg::OperatorNorm<typename ElementTraits<ElementType>::RationalType>(GeNuSys::LinAlg::Algorithms::getJordanForm(invM)))
        {
        }

        template<typename ElementType>
        const GeNuSys::LinAlg::Matrix<ElementType>& RadixProperties<ElementType>::getBase() const
        {
            return M;
        }

        template<typename ElementType>
        const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& RadixProperties<ElementType>::getInverse() const
        {
            return invM;
        }

        template<typename ElementType>
        const GeNuSys::LinAlg::Matrix<ElementType>& RadixProperties<ElementType>::getAdjoint() const
        {
            return adjM;
        }

        template<typename ElementType>
        unsigned int RadixProperties<ElementType>::getSize() const
        {
            return M.getRows();
        }

        template<typename ElementType>
        const ElementType& RadixProperties<ElementType>::getDet() const
        {
            return detM;
        }

        template<typename ElementType>
        const ElementType& RadixProperties<ElementType>::getAbsDet() const
        {
            return absDetM;
        }

        template<typename ElementType>
        const GeNuSys::LinAlg::SmithNormalForm<ElementType>& RadixProperties<ElementType>::getSmithNormalForm() const
        {
            return smithNormalForm;
        }

        template<typename ElementType>
        const GeNuSys::LinAlg::OperatorNorm<typename ElementTraits<ElementType>::RationalType>& RadixProperties<ElementType>::getOperatorNorm() const
        {
            return operatorNorm;
        }

    }
}
