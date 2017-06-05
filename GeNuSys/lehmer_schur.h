/*
GeNuSys - computations with generalized number systems
Copyright (C) 2015-2017  Bence Németh
Copyright (C) 2017  Tamás Krutki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
