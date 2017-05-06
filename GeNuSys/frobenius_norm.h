#ifndef GENUSYS_LINALG_FROBENIUS_NORM_H_
#define GENUSYS_LINALG_FROBENIUS_NORM_H_

#include "element_traits.h"

#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        struct FrobeniusNorm
        {

            template<typename ElementType>
            struct NormType
            {

                typedef typename ElementTraits<typename ElementTraits<ElementType>::AbsSqrType>::RealType Type;

            };

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseMatrix<ElementType>& mat);

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Matrix<ElementType>& mat) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseMatrix<ElementType>& mat) const;

        };

    }
}

// Include implementation
#include "frobenius_norm.hpp"

#endif // GENUSYS_LINALG_FROBENIUS_NORM_H_
