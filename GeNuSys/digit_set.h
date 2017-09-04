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

#ifndef GENUSYS_NUMSYS_DIGIT_SET_H_
#define GENUSYS_NUMSYS_DIGIT_SET_H_

#include <vector>

#include "radix_properties.h"

namespace GeNuSys
{
    namespace NumSys
    {


        struct DigitSet
        {

            template<typename ElementType>
            static std::vector<GeNuSys::LinAlg::SparseVector<ElementType>> getJCanonical(const RadixProperties<ElementType>& props, unsigned int j);

            template<typename ElementType>
            static std::vector<GeNuSys::LinAlg::SparseVector<ElementType>> getJSymmetric(const RadixProperties<ElementType>& props, unsigned int j);

            template<typename ElementType>
            static std::vector<GeNuSys::LinAlg::Vector<ElementType>> getAdjoint(const RadixProperties<ElementType>& props);

            template<typename ElementType>
            static std::vector<GeNuSys::LinAlg::Vector<ElementType>> getDense(const RadixProperties<ElementType>& props);

        };

    }
}

// Include implementation
#include "digit_set.hpp"

#endif // GENUSYS_NUMSYS_DIGIT_SET_H_
