#ifndef GENUSYS_ALGEBRA_LEHMER_SCHUR_H_
#define GENUSYS_ALGEBRA_LEHMER_SCHUR_H_

#include "polynomial.h"

namespace GeNuSys
{
    namespace Algebra
    {

        template<typename ElementType>
        bool lehmerSchur(const Polynomial<ElementType>& poly);

    }
}

// Include implementation
#include "lehmer_schur.hpp"

#endif // GENUSYS_ALGEBRA_LEHMER_SCHUR_H_
