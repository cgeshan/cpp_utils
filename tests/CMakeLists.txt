# ./cpp_utils/test/CMakeLists.txt

include_directories(../src)
add_executable(MatrixTest MatrixTests.cpp)
target_link_libraries(MatrixTest cgClassLib)


add_test(NAME MatrixTest COMMAND cmake --build . --config Release --target MatrixTest)

set_tests_properties(MatrixTest


                    PROPERTIES LABELS "cgClass")