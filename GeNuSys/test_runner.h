#ifndef GENUSYS_TESTS_TEST_RUNNER_H_
#define GENUSYS_TESTS_TEST_RUNNER_H_

#include <string>
#include <iostream>
#include <vector>

#include "test_suite.h"

namespace GeNuSys
{
    namespace Tests
    {

        class TestRunner
        {

            private:

                unsigned int cnt;

                unsigned int passed;

                unsigned int failed;

                std::vector<TestSuite*> suites;

            public:

                TestRunner();

                virtual ~TestRunner();

                void addTestSuite(TestSuite* suite);

                void run();

                int returnCode();

        };

    }
}

#endif // GENUSYS_TESTS_TEST_RUNNER_H_
