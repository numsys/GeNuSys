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
#include <GeNuSys/simultan.h>

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
