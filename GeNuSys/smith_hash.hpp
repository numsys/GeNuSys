#include "element_traits.h"

#include "linalg_traits.h"
#include "linalg_algorithms.h"

namespace GeNuSys
{
    namespace NumSys
    {

        template <
            typename ElementType,
            template<typename> class MatrixType
            >
        SmithHash<ElementType, MatrixType>::SmithHash(const RadixProperties<ElementType>& props)
        {
            const GeNuSys::LinAlg::SmithNormalForm<ElementType>& smithNormalForm = props.getSmithNormalForm();

            unsigned int s = 0;
            while (s < smithNormalForm.S.getRows() && smithNormalForm.S(s, s) == ElementTraits<ElementType>::one())
            {
                ++s;
            }
            size = smithNormalForm.S.getRows() - s;

            U = GeNuSys::LinAlg::Traits::getRows(smithNormalForm.U, s, smithNormalForm.U.getRows());
            G = GeNuSys::LinAlg::Vector<ElementType>(size);
            for (unsigned int i = s; i < smithNormalForm.S.getRows(); ++i)
            {
                G.set(i - s, smithNormalForm.S(i, i));
            }

            prodG = GeNuSys::LinAlg::Vector<ElementType>(size);
            ElementType prod = ElementTraits<ElementType>::one();
            for (unsigned int i = 0; i < size; ++i)
            {
                prodG.set(i, prod);
                prod *= G[i];
            }
        }

        template <
            typename ElementType,
            template<typename> class MatrixType
            >
        GeNuSys::LinAlg::Vector<ElementType> SmithHash<ElementType, MatrixType>::createCache() const
        {
            return GeNuSys::LinAlg::Vector<ElementType>(U.getRows());
        }

        template <
            typename ElementType,
            template<typename> class MatrixType
            >
        ElementType SmithHash<ElementType, MatrixType>::operator()(const GeNuSys::LinAlg::Vector<ElementType>& z) const
        {
            return (*this)(z, createCache());
        }

        template <
            typename ElementType,
            template<typename> class MatrixType
            >
        ElementType SmithHash<ElementType, MatrixType>::operator()(const GeNuSys::LinAlg::Vector<ElementType>& z, GeNuSys::LinAlg::Vector<ElementType>& Uz) const
        {
            GeNuSys::LinAlg::Operations::mat_mul(U, z, Uz);
            ElementType sum = ElementTraits<ElementType>::zero();
            for (unsigned int i = 0; i < size; ++i)
            {
                sum += ElementTraits<ElementType>::mod(Uz[i], G[i]) * prodG[i];
            }

            return sum;
        }

    }
}
