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
