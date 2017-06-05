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

#ifndef GENUSYS_TESTS_TEST_UTILS_H_
#define GENUSYS_TESTS_TEST_UTILS_H_

#include <GeNuSys/element_traits.h>

#include <GeNuSys/vector.h>
#include <GeNuSys/sparse_vector.h>
#include <GeNuSys/matrix.h>
#include <GeNuSys/sparse_matrix.h>

namespace GeNuSys
{
    namespace Tests
    {

        struct TestUtils
        {

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB);

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB);

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB);

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB);

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB);

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB);

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB);

            template<typename T, typename S>
            static bool equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB);

        };

    }
}

// Include implementation
#include "test_utils.hpp"

#endif // GENUSYS_TESTS_TEST_UTILS_H_
