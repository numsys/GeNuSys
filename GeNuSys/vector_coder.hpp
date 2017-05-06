#include "element_traits.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template<typename ElementType>
        unsigned long long VectorCoder::encode(const GeNuSys::LinAlg::Vector<ElementType>& z, bool& valid)
        {
            valid = true;
            unsigned long long code = 0;
            for (int j = z.getLength() - 1; j > 0; --j)
            {
                if (lowerBound[j] > ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[j])
                        || ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[j]) > upperBound[j])
                {
                    valid = false;
                }
                code += ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[j]) - lowerBound[j];
                code *= varBase[j - 1];
            }
            code += ElementTraits<ElementType>::template asTypeUnsafe<long long>(z[0]) - lowerBound[0];

            return code;
        }

        template<typename ElementType>
        void VectorCoder::decode(unsigned long long code, GeNuSys::LinAlg::Vector<ElementType>& z)
        {
            for (unsigned int j = 0; j < z.getLength(); ++j)
            {
                z.set(j, ElementTraits<long long>::template asTypeUnsafe<ElementType>(code % varBase[j] + lowerBound[j]));
                code /= varBase[j];
            }
        }

    }
}
