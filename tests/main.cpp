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

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include <GeNuSys/element_traits.h>

#include <GeNuSys/vector.h>
#include <GeNuSys/sparse_vector.h>
#include <GeNuSys/matrix.h>
#include <GeNuSys/sparse_matrix.h>

#include <GeNuSys/p_norm.h>
#include <GeNuSys/frobenius_norm.h>
#include <GeNuSys/operator_norm.h>

#include <GeNuSys/linalg_algorithms.h>

#include <GeNuSys/polynomial.h>
#include <GeNuSys/lehmer_schur.h>

#include "test_suite.h"
#include "test_runner.h"
#include "tests.h"

#include <GeNuSys/number_system.h>
#include <GeNuSys/radix_properties.h>
#include <GeNuSys/digit_set.h>
#include <GeNuSys/smith_hash.h>
#include <GeNuSys/numsys_traits.h>
#include <GeNuSys/simultaneous.h>

int runTests()
{
    GeNuSys::Tests::TestRunner testRunner;
    testRunner.addTestSuite(new VectorTest());
    testRunner.addTestSuite(new SparseVectorTest());
    testRunner.addTestSuite(new VectorCmpTest());
    testRunner.addTestSuite(new SparseVectorCmpTest());
    testRunner.addTestSuite(new VectorOperatorTest());
    testRunner.addTestSuite(new SparseVectorOperatorTest());
    testRunner.addTestSuite(new VectorNormTest());
    testRunner.addTestSuite(new MatrixTest());
    testRunner.addTestSuite(new MatrixNormTest());
    testRunner.run();

    return testRunner.returnCode();
}

int main()
{
    int testsResult = runTests();
    return testsResult;
}
