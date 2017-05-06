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
