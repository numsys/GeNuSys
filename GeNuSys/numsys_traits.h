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

#ifndef GENUSYS_NUMSYS_TRAITS_H_
#define GENUSYS_NUMSYS_TRAITS_H_

#include "element_traits.h"

#include "matrix.h"
#include "vector.h"

#include "digit_set.h"

namespace GeNuSys
{
    namespace NumSys
    {

        struct Traits
        {

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static void getBounds(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                                  const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& T,
                                  std::vector<int>& lowerBound, std::vector<int>& upperBound);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static void getBounds(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                                  std::vector<int>& lowerBound, std::vector<int>& upperBound);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static unsigned long long getVolume(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static unsigned long long getVolume(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                                                const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& T);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static GeNuSys::LinAlg::Matrix<ElementType> findBasisTransformation(
                const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                const unsigned int candNum, const unsigned int mutateNum, const unsigned int noImprLimit);

        };

    }
}

// Include implementation
#include "numsys_traits.hpp"

#endif // GENUSYS_NUMSYS_TRAITS_H_
