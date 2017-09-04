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

#ifndef GENUSYS_NUMSYS_BASE_PROPERTIES_H_
#define GENUSYS_NUMSYS_BASE_PROPERTIES_H_

#include "matrix.h"
#include "linalg_algorithms.h"
#include "operator_norm.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template<typename ElementType>
        class RadixProperties
        {

            private:

                GeNuSys::LinAlg::Matrix<ElementType> M;

                GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType> invM;

                GeNuSys::LinAlg::Matrix<ElementType> adjM;

                ElementType detM;

                ElementType absDetM;

                GeNuSys::LinAlg::SmithNormalForm<ElementType> smithNormalForm;

                GeNuSys::LinAlg::OperatorNorm<typename ElementTraits<ElementType>::RationalType> operatorNorm;

            public:

                RadixProperties(const GeNuSys::LinAlg::Matrix<ElementType>& M);

                const GeNuSys::LinAlg::Matrix<ElementType>& getBase() const;

                const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& getInverse() const;

                const GeNuSys::LinAlg::Matrix<ElementType>& getAdjoint() const;

                unsigned int getSize() const;

                const ElementType& getDet() const;

                const ElementType& getAbsDet() const;

                const GeNuSys::LinAlg::SmithNormalForm<ElementType>& getSmithNormalForm() const;

                const GeNuSys::LinAlg::OperatorNorm<typename ElementTraits<ElementType>::RationalType>& getOperatorNorm() const;

        };

    }
}

// Include implementation
#include "radix_properties.hpp"

#endif // GENUSYS_NUMSYS_BASE_PROPERTIES_H_
