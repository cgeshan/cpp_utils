# ./cpp_utils/src/CMakeLists.txt

file(GLOB SRC "*.cpp")
file(GLOB HDRS "*.h")

add_library(cgClassLib ${SRC} ${HDRS})
