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
