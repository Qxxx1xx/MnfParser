#include <test_siso_beam.h>
#include <test_siso_block.h>
#include <test_two_input_single_output_block_beam.h>

// main 函数，设置并运行所有测试
int
main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
