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
