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
        unsigned long long VectorCoder::encode(const GeNuSys::LinAlg::Vector<ElementType>& z, bool& valid)
        {
            valid = true;
            unsigned long long code = 0;
            for (int j = z.getLength() - 1; j > 0; --j)
            {
                if ((lowerBound[j] > ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[j]))
                        || (ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[j]) > upperBound[j]))
                {
                    valid = false;
                }
                code += ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[j]) - lowerBound[j];
                code *= varBase[j - 1];
            }
            if ((lowerBound[0] > ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[0]))
                    || (ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[0]) > upperBound[0]))
            {
                valid = false;
            }
            code += ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[0]) - lowerBound[0];

            return code;
        }

        template<typename ElementType>
        void VectorCoder::decode(unsigned long long code, GeNuSys::LinAlg::Vector<ElementType>& z)
        {
            for (unsigned int j = 0; j < z.getLength(); ++j)
            {
                z.set(j, ElementTraits<long long>::template asTypeUnsafe<ElementType>(code % varBase[j] + lowerBound[j]));
                code /= varBase[j];
            }
        }

    }
}
