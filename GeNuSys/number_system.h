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

#ifndef GENUSYS_NUMSYS_NUMBER_SYSTEM_H_
#define GENUSYS_NUMSYS_NUMBER_SYSTEM_H_

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

#include "linalg_algorithms.h"

#include "radix_properties.h"
#include "digit_set.h"
#include "hash_table.h"
#include "smith_hash.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType,
            typename Norm
            >
        class NumberSystem
        {

            private:

                RadixProperties<ElementType> props;

                std::vector<VectorType<ElementType>> digitSet;

                SmithHash<ElementType, MatrixType> hash;

                HashTable<ElementType, VectorType, MatrixType> hashTable;

                Norm norm;

            public:

                NumberSystem(const RadixProperties<ElementType>& props, const std::vector<VectorType<ElementType>>& digitSet, const Norm& norm);

                GeNuSys::LinAlg::Vector<ElementType> phi(const GeNuSys::LinAlg::Vector<ElementType>& z) const;

                const VectorType<ElementType>& phi(GeNuSys::LinAlg::Vector<ElementType>& z, GeNuSys::LinAlg::Vector<ElementType>& phiZ, GeNuSys::LinAlg::Vector<ElementType>& Uz) const;

                std::vector<VectorType<ElementType>> getExpansion(const GeNuSys::LinAlg::Vector<ElementType>& z);

                std::vector<GeNuSys::LinAlg::Vector<ElementType>> getOrbit(const GeNuSys::LinAlg::Vector<ElementType>& z);

                std::vector<std::vector<GeNuSys::LinAlg::Vector<ElementType>>> getCycles();

        };

    }
}

// Include implementation
#include "number_system.hpp"

#endif
