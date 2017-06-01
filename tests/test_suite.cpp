#include "test_suite.h"

namespace GeNuSys
{
    namespace Tests
    {

        TestSuite::TestSuite(const std::string& name): name(name), cnt(0), passed(0), failed(0)
        {
        }

        TestSuite::~TestSuite()
        {
        }

        void TestSuite::fail(const std::string& message)
        {
            ++cnt;
            ++failed;
            std::cout << "ASSERT " << cnt << ": " << (message.length() == 0 ? "" : message + " - ") << "FAILED" << std::endl;
        }

        void TestSuite::pass(const std::string& message)
        {
            ++cnt;
            ++passed;
            std::cout << "ASSERT " << cnt << ": " << (message.length() == 0 ? "" : message + " - ") << "PASSED" << std::endl;
        }

        void TestSuite::runTestSuite()
        {
            std::cout << "============================================================" << std::endl;
            std::cout << "RUNNING TEST SUITE - " << name << std::endl;
            std::cout << "===========================================================" << std::endl;
            run();
            std::cout << "===========================================================" << std::endl;
            if (failed == 0)
            {
                std::cout << "RUN " << cnt << " TEST(S) - ALL TESTS PASSED!" << std::endl;
            }
            else
            {
                std::cout << "RUN " << cnt << " TEST(S) - FAILED: " << failed << std::endl;
            }
            std::cout << std::endl;
        }

        void TestSuite::assertTrue(bool value, const std::string& message)
        {
            if (value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: TRUE - Actual value: FALSE" << std::endl;
            }
        }

        void TestSuite::assertFalse(bool value, const std::string& message)
        {
            if (!value)
            {
                pass(message);
            }
            else
            {
                fail(message);
                std::cout << "Expected: FALSE - Actual value: TRUE" << std::endl;
            }
        }

    }
}
