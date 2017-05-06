#ifndef GENUSYS_NUMSYS_TRAITS_H_
#define GENUSYS_NUMSYS_TRAITS_H_

#include "element_traits.h"

#include "matrix.h"
#include "vector.h"

#include "digit_set.h"

namespace GeNuSys
{
    namespace NumSys
    {

        struct Traits
        {

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static void getBounds(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                                  const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& T,
                                  std::vector<int>& lowerBound, std::vector<int>& upperBound);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static void getBounds(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                                  std::vector<int>& lowerBound, std::vector<int>& upperBound);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static unsigned long long getVolume(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static unsigned long long getVolume(const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                                                const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& T);

            template <
                typename ElementType,
                template<typename> class VectorType
                >
            static GeNuSys::LinAlg::Matrix<ElementType> findBasisTransformation(
                const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& invM, const std::vector<VectorType<ElementType>>& digitSet,
                const unsigned int candNum, const unsigned int mutateNum, const unsigned int noImprLimit);

        };

    }
}

// Include implementation
#include "numsys_traits.hpp"

#endif // GENUSYS_NUMSYS_TRAITS_H_
