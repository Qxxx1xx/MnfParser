# 查找 GoogleTest
find_package(GTest REQUIRED)

# 给MnfParserLib添加测试
file(GLOB_RECURSE srcs CONFIGURE_DEPENDS *.cpp)
file(GLOB_RECURSE incs CONFIGURE_DEPENDS *.h)

add_executable(MnfParserTests ${srcs} ${incs})
target_compile_definitions(MnfParserTests PRIVATE "MnfParserTESTS_DATA_DIR=\"${CMAKE_SOURCE_DIR}/tests/data/\"")
target_include_directories(MnfParserTests PRIVATE .)
target_link_libraries(MnfParserTests PUBLIC MnfParserLib GTest::gtest GTest::gtest_main)
add_test(NAME MnfParserTests COMMAND MnfParserTests)
