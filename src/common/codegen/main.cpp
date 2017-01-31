#include <string>
#include "outputs/outputs.hpp"

#ifdef SP2_ENABLE_TESTS
#include <gtest/gtest.h>
#endif // SP2_ENABLE_TESTS

int run_codegen(int argc, char *argv[]);

int main(int argc, char *argv[])
{
#ifdef SP2_ENABLE_TESTS
    for (int i = 0; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--test")
        {
            ::testing::InitGoogleTest(&argc, argv);
            return RUN_ALL_TESTS();
        }
    }
#endif // SP2_ENABLE_TESTS

    // generate exp_table.hpp
    exp_table_hpp();
}