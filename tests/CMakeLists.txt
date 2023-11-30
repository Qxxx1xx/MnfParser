# 查找 GoogleTest
find_package(GTest REQUIRED)

# 给MnfParserLib添加测试
add_executable(MnfParserTests test_MnfParserLib.cpp)
target_compile_definitions(MnfParserTests PRIVATE "MnfParserTESTS_DATA_DIR=\"${CMAKE_SOURCE_DIR}/tests/data/\"")
target_link_libraries(MnfParserTests PUBLIC MnfParserLib GTest::gtest GTest::gtest_main)
add_test(NAME MnfParserTests COMMAND MnfParserTests)
