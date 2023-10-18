# ./cpp_utils/test/CMakeLists.txt

include_directories(../src)

add_executable(MatrixTest MatrixTests.cpp)
target_link_libraries(MatrixTest cgClassLib)

add_executable(VectorTest VectorTests.cpp)
target_link_libraries(VectorTest cgClassLib)




add_test(NAME MatrixTest COMMAND cmake --build . --config Release --target MatrixTest)
add_test(NAME VectorTest COMMAND cmake --build . --config Release --target VectorTest)

set_tests_properties(MatrixTest
                     VectorTest

                    PROPERTIES LABELS "cgClass")