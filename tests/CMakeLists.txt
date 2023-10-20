# ./cpp_utils/test/CMakeLists.txt

include_directories(../src)

add_executable(MatrixTest MatrixTests.cpp)
target_link_libraries(MatrixTest cgClassLib)

add_executable(VectorTest VectorTests.cpp)
target_link_libraries(VectorTest cgClassLib)




add_test(NAME BuildAll COMMAND cmake --build .)
add_test(NAME RunMatrixTest COMMAND ./MatrixTest)
add_test(NAME RunVectorTest COMMAND ./VectorTest)

set_tests_properties(BuildAll
                    
                    RunMatrixTest
                    RunVectorTest

                    PROPERTIES LABELS "cgClass")