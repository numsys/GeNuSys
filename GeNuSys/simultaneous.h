#ifndef GENUSYS_NUMSYS_SIMULTANEOUS_H_
#define GENUSYS_NUMSYS_SIMULTANEOUS_H_

#include <vector>

#include "element_traits.h"

#include "vector.h"
#include "p_norm.h"

#include "radix_properties.h"
#include "smith_hash.h"

#include <fstream>

namespace GeNuSys
{
    namespace NumSys
    {

        struct Simultaneous
        {

            template<typename ElementType>
            static GeNuSys::LinAlg::Matrix<ElementType> getSimultaneous(ElementType a, ElementType b);

            template<typename ElementType, typename Norm>
            static void tryDigitCandidate(const ElementType i, const ElementType j, GeNuSys::LinAlg::Vector<ElementType>& vct,
                                          std::vector<bool>& digitFound, std::vector<typename Norm::template NormType<ElementType>::Type>& digitNorm, std::vector<GeNuSys::LinAlg::Vector<ElementType>>& digits,
                                          int& numFound,
                                          const Norm& norm, typename Norm::template NormType<ElementType>::Type& maxNorm, bool& flag,
                                          const SmithHash<ElementType, GeNuSys::LinAlg::Matrix>& hash, GeNuSys::LinAlg::Vector<ElementType>& Uz);

            template<typename ElementType, typename Norm>
            static std::vector<GeNuSys::LinAlg::Vector<ElementType>> getDigitSet(const RadixProperties<ElementType>& props, const Norm& norm, bool global);

            /*
            template<typename ElementType>
            static std::vector<GeNuSys::LinAlg::Vector<ElementType>> readDigitSet()
            {
                std::vector<GeNuSys::LinAlg::Vector<ElementType>> digits;
                ElementType a[147] = { 4, -4, -4, 6, -3, 0, 1, 2, -2, -2, -2, 1, -3, 3, 0, 6, 2, 2, -4, -1, -1, -2, 1, 3, -3, 7, 6, -5, 3, -4, -8, 3, -1, -1, 8, 4, 1, 3, 1, -3, -7, 7, 0, 2, -2, -1, -5, 0, 0, -3, -2, -6, 1, -1, -3, 7, 2, -7, 2, -2, 2, 5, -4, 1, -6, -4, -2, -6, 1, 4, -4, -1, -6, 4, -3, 4, 5, 0, 1, 0, 0, -2, 3, -1, 2, 2, 5, 6, 1, 4, -2, 1, 1, -3, -3, 5, 0, 4, 0, 5, 1, -4, -5, -1, 4, -4, -1, -3, -5, 5, 6, 0, -5, -1, 0, -6, -3, 1, 0, -4, -8, -1, 4, 2, -2, 0, 3, -5, 2, -2, 2, 3, 3, -1, 2, -2, 3, -1, -5, -2, -7, -7, 3, 1, -3, 2, 5};
                ElementType b[147] = { 6, -6, -1, 1, 2, -3, 0, 3, 5, 4, 2, -5, -3, 6, 6, 4, 6, -1, -5, -6, -3, 3, -2, -2, -7, 2, 2, -4, 0, 1, -1, 1, 3, -2, 2, -1, 5, 5, -1, -2, 0, 1, 3, -3, -8, -5, -3, -5, 5, 0, 1, -1, 6, 0, -6, 3, 7, -1, -5, -3, 0, 5, 2, 1, 0, 0, -5, -3, 7, 1, 3, -4, -4, 5, -5, -2, 2, 7, -4, 0, -2, -4, -3, -1, 2, 8, 3, 3, 8, 3, 0, -3, 2, 4, -1, 0, -1, 0, 2, 1, 3, -3, 1, 5, 4, -2, -7, 1, -5, 4, 0, -4, -2, 2, 1, -2, 3, 4, -6, -4, -2, 1, 2, -4, -2, 4, -1, 0, -2, -7, 1, 3, 2, 4, 5, -1, 4, 6, -1, -6, -3, -2, 7, -6, -4, 4, -1 };
                GeNuSys::LinAlg::Vector<ElementType> digit(6);
                for (int i = 0; i < 147; ++i)
                {
                    digit.set(0, a[i]);
                    digit.set(1, b[i]);
                    digit.set(2, a[i]);
                    digit.set(3, b[i]);
                    digit.set(4, a[i]);
                    digit.set(5, b[i]);
                    digits.push_back(digit);
                }
                return digits;
            }
            */

        };

    }
}

#include "simultaneous.hpp"

#endif // GENUSYS_NUMSYS_SIMULTANEOUS_H_
