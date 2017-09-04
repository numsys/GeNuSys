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

#ifndef GENUSYS_NUMSYS_HASH_TABLE_H_
#define GENUSYS_NUMSYS_HASH_TABLE_H_

#include <vector>

#include "smith_hash.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType
            >
        class HashTable
        {

            private:

                SmithHash<ElementType, MatrixType> hash;

                std::vector<VectorType<ElementType>> digitSet;

                std::vector<VectorType<ElementType>> hashTable;

            public:

                HashTable(const SmithHash<ElementType, MatrixType>& hash, const std::vector<VectorType<ElementType>>& digitSet);

                const VectorType<ElementType>& operator()(const GeNuSys::LinAlg::Vector<ElementType>& z) const;

                const VectorType<ElementType>& operator()(const GeNuSys::LinAlg::Vector<ElementType>& z, GeNuSys::LinAlg::Vector<ElementType>& Uz) const;

        };

    }
}

// Include implementation
#include "hash_table.hpp"

#endif // GENUSYS_NUMSYS_HASH_TABLE_H_
