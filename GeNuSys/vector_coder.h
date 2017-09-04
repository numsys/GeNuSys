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

#ifndef GENUSYS_NUMSYS_VECTOR_CODER_H_
#define GENUSYS_NUMSYS_VECTOR_CODER_H_

#include <vector>

#include "vector.h"

namespace GeNuSys
{
    namespace NumSys
    {

        class VectorCoder
        {

            private:

                std::vector<int> lowerBound;

                std::vector<int> upperBound;

                std::vector<unsigned int> varBase;

                unsigned long long size;

            public:

                VectorCoder(const std::vector<int>& lowerBound, const std::vector<int>& upperBound): lowerBound(lowerBound), upperBound(upperBound)
                {
                    varBase = std::vector<unsigned int>(lowerBound.size());
                    size = 1;
                    for (unsigned int i = 0; i < lowerBound.size(); ++i)
                    {
                        varBase[i] = upperBound[i] - lowerBound[i] + 1;
                        size *= varBase[i];
                    }
                }

                unsigned long long getSize()
                {
                    return size;
                }

                template<typename ElementType>
                unsigned long long encode(const GeNuSys::LinAlg::Vector<ElementType>& z, bool& valid);

                template<typename ElementType>
                void decode(unsigned long long code, GeNuSys::LinAlg::Vector<ElementType>& z);

        };

    }
}

// Include implementation
#include "vector_coder.hpp"

#endif // GENUSYS_NUMSYS_VECTOR_CODER_H_
