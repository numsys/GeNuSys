#ifndef GENUSYS_LINALG_P_NORM_H_
#define GENUSYS_LINALG_P_NORM_H_

#include "element_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys
{
    namespace LinAlg
    {

        template<unsigned int P>
        struct PNorm
        {

            template<typename ElementType>
            struct NormType
            {

                typedef typename ElementTraits<typename ElementTraits<ElementType>::AbsType>::RealType Type;

            };

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Vector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseVector<ElementType>& vct);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const Matrix<ElementType>& mat);

            template<typename ElementType>
            static typename NormType<ElementType>::Type norm(const SparseMatrix<ElementType>& mat);

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Vector<ElementType>& vct) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseVector<ElementType>& vct) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const Matrix<ElementType>& mat) const;

            template<typename ElementType>
            typename NormType<ElementType>::Type operator()(const SparseMatrix<ElementType>& mat) const;

        };

    }
}

// Include implementation
#include "p_norm.hpp"

#endif // GENUSYS_LINALG_P_NORM_H_
