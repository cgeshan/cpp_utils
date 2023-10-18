# ./cpp_utils/CMakeLists.txt

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
cmake_minimum_required(VERSION 3.15)

project(cpp_utils)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

