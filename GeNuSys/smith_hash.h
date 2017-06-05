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

#ifndef GENUSYS_NUMSYS_SMITH_HASH_H_
#define GENUSYS_NUMSYS_SMITH_HASH_H_

#include "vector.h"

#include "vector.h"
#include "radix_properties.h"
#include "digit_set.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class MatrixType
            >
        class SmithHash
        {

            private:

                unsigned int size;

                MatrixType<ElementType> U;

                GeNuSys::LinAlg::Vector<ElementType> G;

                GeNuSys::LinAlg::Vector<ElementType> prodG;

            public:

                SmithHash(const RadixProperties<ElementType>& props);

                GeNuSys::LinAlg::Vector<ElementType> createCache() const;

                ElementType operator()(const GeNuSys::LinAlg::Vector<ElementType>& z) const;

                ElementType operator()(const GeNuSys::LinAlg::Vector<ElementType>& z, GeNuSys::LinAlg::Vector<ElementType>& Uz) const;

        };

    }
}

#include "smith_hash.hpp"

#endif // GENUSYS_NUMSYS_SMITH_HASH_H_
