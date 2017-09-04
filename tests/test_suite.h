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

#ifndef GENUSYS_TESTS_TEST_SUITE_H_
#define GENUSYS_TESTS_TEST_SUITE_H_

#include <string>
#include <iostream>

namespace GeNuSys
{
    namespace Tests
    {

        class TestSuite
        {

                friend class TestRunner;

            private:

                std::string name;

                int cnt;

                int passed;

                int failed;

                void fail(const std::string& message);

                void pass(const std::string& message);

            public:

                TestSuite(const std::string& name);

                virtual ~TestSuite();

                void runTestSuite();

                virtual void run() = 0;

                void assertTrue(bool value, const std::string& message = "");

                void assertFalse(bool value, const std::string& message = "");

                template<typename T>
                void assertEqual(const T& expected, const T& value, const std::string& message = "");

                template<typename T>
                void assertNotEqual(const T& expected, const T& value, const std::string& message = "");

                template<typename T>
                void assertLess(const T& expected, const T& value, const std::string& message = "");

                template<typename T>
                void assertLessOrEqual(const T& expected, const T& value, const std::string& message = "");

                template<typename T>
                void assertMore(const T& expected, const T& value, const std::string& message = "");

                template<typename T>
                void assertMoreOrEqual(const T& expected, const T& value, const std::string& message = "");

                template<typename T>
                void assertBetween(const T& low, const T& high, const T& value, const std::string& message = "");

                template<typename T>
                void assertNotBetween(const T& low, const T& high, const T& value, const std::string& message = "");

                template<typename T>
                void assertDifference(const T& expected, const T& threshold, const T& value, const std::string& message = "");

        };

        template<typename T>
        void TestSuite::assertEqual(const T& expected, const T& value, const std::string& message)
        {
            if (expected == value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: " << expected << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertNotEqual(const T& expected, const T& value, const std::string& message)
        {
            if (expected != value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: NOT " << expected << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertLess(const T& expected, const T& value, const std::string& message)
        {
            if (expected > value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: LESS THAN " << expected << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertLessOrEqual(const T& expected, const T& value, const std::string& message)
        {
            if (expected >= value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: LESS THAN OR EQUAL TO " << expected << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertMore(const T& expected, const T& value, const std::string& message)
        {
            if (expected < value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: MORE THAN " << expected << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertMoreOrEqual(const T& expected, const T& value, const std::string& message)
        {
            if (expected <= value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: MORE THAN OR EQUAL TO " << expected << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertBetween(const T& low, const T& high, const T& value, const std::string& message)
        {
            if (low <= value && value <= high)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: BETWEEN " << low << " AND " << high << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertNotBetween(const T& low, const T& high, const T& value, const std::string& message)
        {
            if (!(low <= value && value <= high))
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: NOT BETWEEN " << low << " AND " << high << " - Actual value: " << value << std::endl;
            }
        }

        template<typename T>
        void TestSuite::assertDifference(const T& expected, const T& threshold, const T& value, const std::string& message)
        {
            T diff = expected - value;
            if (-threshold <= diff && diff <= threshold)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: CLOSER TO " << expected << " THAN " << threshold << " - Actual value: " << value << " (DIFFERENCE: " << diff << ")" << std::endl;
            }
        }

    }
}

#endif // GENUSYS_TESTS_TEST_SUITE_H_
