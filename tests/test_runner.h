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
