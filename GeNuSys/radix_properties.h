#ifndef GENUSYS_NUMSYS_BASE_PROPERTIES_H_
#define GENUSYS_NUMSYS_BASE_PROPERTIES_H_

#include "matrix.h"
#include "linalg_algorithms.h"
#include "operator_norm.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template<typename ElementType>
        class RadixProperties
        {

            private:

                GeNuSys::LinAlg::Matrix<ElementType> M;

                GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType> invM;

                GeNuSys::LinAlg::Matrix<ElementType> adjM;

                ElementType detM;

                ElementType absDetM;

                GeNuSys::LinAlg::SmithNormalForm<ElementType> smithNormalForm;

                GeNuSys::LinAlg::OperatorNorm<typename ElementTraits<ElementType>::RationalType> operatorNorm;

            public:

                RadixProperties(const GeNuSys::LinAlg::Matrix<ElementType>& M);

                const GeNuSys::LinAlg::Matrix<ElementType>& getBase() const;

                const GeNuSys::LinAlg::Matrix<typename ElementTraits<ElementType>::RationalType>& getInverse() const;

                const GeNuSys::LinAlg::Matrix<ElementType>& getAdjoint() const;

                unsigned int getSize() const;

                const ElementType& getDet() const;

                const ElementType& getAbsDet() const;

                const GeNuSys::LinAlg::SmithNormalForm<ElementType>& getSmithNormalForm() const;

                const GeNuSys::LinAlg::OperatorNorm<typename ElementTraits<ElementType>::RationalType>& getOperatorNorm() const;

        };

    }
}

// Include implementation
#include "radix_properties.hpp"

#endif // GENUSYS_NUMSYS_BASE_PROPERTIES_H_
