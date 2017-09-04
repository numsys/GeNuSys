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

#include "test_runner.h"

namespace GeNuSys
{
    namespace Tests
    {

        TestRunner::TestRunner(): cnt(0), passed(0), failed(0)
        {
        }

        TestRunner::~TestRunner()
        {
            for (unsigned int i = 0; i < suites.size(); ++i)
            {
                delete suites[i];
            }
        }

        void TestRunner::addTestSuite(TestSuite* suite)
        {
            suites.push_back(suite);
        }

        void TestRunner::run()
        {
            for (unsigned int i = 0; i < suites.size(); ++i)
            {
                suites[i]->runTestSuite();
                cnt += suites[i]->cnt;
                passed += suites[i]->passed;
                failed += suites[i]->failed;
            }

            std::cout << std::endl;
            std::cout << "===========================================================" << std::endl;
            for (unsigned int i = 0; i < suites.size(); ++i)
            {
                std::cout << "TEST SUITE " << suites[i]->name << " | " << suites[i]->cnt << " | " << suites[i]->passed << " | " << suites[i]->failed << std::endl;
            }
            std::cout << "===========================================================" << std::endl;
            if (failed == 0)
            {
                std::cout << "RUN " << cnt << " TEST(S) IN " << suites.size() << " SUITES - ALL TESTS PASSED!" << std::endl;
            }
            else
            {
                std::cout << "RUN " << cnt << " TEST(S) IN " << suites.size() << " SUITES - FAILED: " << failed << std::endl;
            }
            std::cout << std::endl;
            std::cout << "===========================================================" << std::endl;
        }

        int TestRunner::returnCode()
        {
            return (failed == 0) ? 0 : 1;
        }

    }
}
