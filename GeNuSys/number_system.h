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
