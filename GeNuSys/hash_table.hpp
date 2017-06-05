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

namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType
            >
        HashTable<ElementType, VectorType, MatrixType>::HashTable(const SmithHash<ElementType, MatrixType>& hash, const std::vector<VectorType<ElementType>>& digitSet): hash(hash), digitSet(digitSet)
        {
            // TODO: VERIFY CRS PROPERTY!!!
            hashTable = std::vector<VectorType<ElementType>>(digitSet.size());
            GeNuSys::LinAlg::Vector<ElementType> Uz = hash.createCache();
            for (unsigned int i = 0; i < digitSet.size(); ++i)
            {
                hashTable[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(hash(digitSet[i], Uz))] = digitSet[i];
            }
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType
            >
        const VectorType<ElementType>& HashTable<ElementType, VectorType, MatrixType>::operator()(const GeNuSys::LinAlg::Vector<ElementType>& z) const
        {
            return hashTable[hash(z, hash.createCache())];
        }

        template <
            typename ElementType,
            template<typename> class VectorType,
            template<typename> class MatrixType
            >
        const VectorType<ElementType>& HashTable<ElementType, VectorType, MatrixType>::operator()(const GeNuSys::LinAlg::Vector<ElementType>& z, GeNuSys::LinAlg::Vector<ElementType>& Uz) const
        {
            return hashTable[ElementTraits<ElementType>::template asTypeUnsafe<unsigned long int>(hash(z, Uz))];
        }

    }
}
