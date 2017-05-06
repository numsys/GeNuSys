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
