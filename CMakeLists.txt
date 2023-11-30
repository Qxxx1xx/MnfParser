# 指定最低所需的 CMake 版本
cmake_minimum_required(VERSION 3.15)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# 指定C++标准为C++20
set(CMAKE_CXX_STANDARD 20)

# 如果编译器不支持C++17将会报错
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 不使用任何编译器特定的扩展
set(CMAKE_CXX_EXTENSIONS OFF)

# 执行安装vcpkg
execute_process(
    COMMAND cmd /c "${CMAKE_CURRENT_SOURCE_DIR}/install_vcpkg.bat"
)
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake"
    CACHE STRING "Vcpkg toolchain file")
project(MnfParser)

# 启用测试支持
enable_testing()

find_package(Eigen3 CONFIG REQUIRED)

# 将子项目MnfParser添加进来
add_subdirectory(MnfParser)
add_subdirectory(tests)