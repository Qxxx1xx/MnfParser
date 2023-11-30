# 定义库源文件
file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)

# 定义主程序源文件
set(main_src src/main.cpp)

# 排除主程序文件从库源文件
list(REMOVE_ITEM srcs ${main_src})

# 添加静态库目标
add_library(MnfParserLib STATIC ${srcs})
target_include_directories(MnfParserLib PUBLIC include)
target_link_libraries(MnfParserLib PUBLIC Eigen3::Eigen)

# 添加可执行文件目标，并链接静态库
add_executable(MnfParser ${main_src})
target_link_libraries(MnfParser PUBLIC MnfParserLib)
