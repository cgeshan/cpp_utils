#include <iostream>
#include "../src/Matrix.h"

int main()
{
    // Test 1: Matrix construction and element access
    std::cout << "Test 1: Matrix construction and element access..." << std::endl;
    Matrix<int> mat(2, 3);
    assert(mat.getRowSize() == 2);
    assert(mat.getColSize() == 3);

    Matrix<double> squareMatrix(8);
    assert(squareMatrix.getRowSize() == 8);
    assert(squareMatrix.getColSize() == 8);

    // Test 2: Copy constructor
    std::cout << "Test 2: Copy constructor..." << std::endl;
    Matrix<int> mat2(mat);
    assert(mat2.getRowSize() == 2);
    assert(mat2.getColSize() == 3);

    // Test 3: Matrix construction with vector of values
    std::cout << "Test 3: Matrix construction with vector of values..." << std::endl;
    std::vector<std::vector<int>> values = {{1, 2, 3}, {4, 5, 6}};
    Matrix<int> mat3(values);
    assert(mat3.getRowSize() == 2);
    assert(mat3.getColSize() == 3);

    // Test 4: Matrix construction with specified dimensions and vector of values
    std::cout << "Test 4: Matrix construction with specified dimensions and vector of values..." << std::endl;
    std::vector<int> values2 = {1, 2, 3, 4, 5, 6};
    Matrix<int> mat4(2, 3, values2);
    assert(mat4.getRowSize() == 2);
    assert(mat4.getColSize() == 3);
    assert(mat4(0, 1) == 2);
    assert(mat4(1, 2) == 6);

    // Test 5: Copy assignment operator
    std::cout << "Test 5: Copy assignment operator..." << std::endl;
    Matrix<int> mat5 = mat4;
    assert(mat5 == mat4);

    // Test 6: Equality operator
    std::cout << "Test 6: Equality operator..." << std::endl;
    assert(mat4 == mat5);

    // Test 7: Addition operator
    std::cout << "Test 7: Addition operator..." << std::endl;
    Matrix<int> mat6({{1, 2}, {3, 4}});
    Matrix<int> mat7({{5, 6}, {7, 8}});
    Matrix<int> sum = mat6 + mat7;
    assert(sum(0, 1) == 8);
    assert(sum(1, 1) == 12);

    // Test 8: Subtraction operator
    std::cout << "Test 8: Subtraction operator..." << std::endl;
    Matrix<int> diff = mat6 - mat7;
    assert(diff(0, 1) == -4);
    assert(diff(1, 1) == -4);

    // Test 9: Multiplication operator
    std::cout << "Test 9: Multiplication operator..." << std::endl;
    Matrix<int> product = mat6 * mat7;
    assert(product(0, 0) == 19);
    assert(product(1, 1) == 50);

    // Test 10: Multiplication by scalar
    std::cout << "Test 10: Multiplication by scalar..." << std::endl;
    Matrix<int> scaled = mat4 * 2;
    assert(scaled(0, 1) == 4);
    assert(scaled(1, 2) == 12);

    // Test 11: Division by scalar
    std::cout << "Test 11: Division by scalar..." << std::endl;
    Matrix<int> divided = mat4 / 2;
    assert(divided(0, 1) == 1);
    assert(divided(1, 2) == 3);

    // Test 12: Addition assignment operator
    std::cout << "Test 12: Addition assignment operator..." << std::endl;
    Matrix<int> mat8(mat4);
    mat8 += mat5;
    assert(mat8(0, 1) == 4);
    assert(mat8(1, 2) == 12);

    // Test 13: Subtraction assignment operator
    std::cout << "Test 13: Subtraction assignment operator..." << std::endl;
    Matrix<int> mat9(mat4);
    mat9 -= mat5;
    assert(mat9(0, 1) == 0);
    assert(mat9(1, 2) == 0);

    // Test 14: Multiplication assignment operator
    std::cout << "Test 14: Multiplication assignment operator..." << std::endl;
    Matrix<int> mat10({{1, 2}, {3, 4}});
    Matrix<int> mat11({{5, 6}, {7, 8}});
    mat10 *= mat11;
    assert(mat10 == Matrix<int>({{19, 22}, {43, 50}}));

    // Test 15: Multiplication assignment operator for matrices of different sizes
    std::cout << "Test 15: Multiplication assignment operator for matrices of different sizes..." << std::endl;
    Matrix<int> mat12({{1, 2, 3}, {4, 5, 6}});
    Matrix<int> mat13({{7, 8}, {9, 10}, {11, 12}});
    mat12 *= mat13;
    assert(mat12 == Matrix<int>({{58, 64}, {139, 154}}));

    // Test 16: Scalar multiplication assignment operator
    std::cout << "Test 16: Scalar multiplication assignment operator..." << std::endl;

    // Test 16.1: Scalar multiplication for a 2x2 matrix
    Matrix<int> mat14({{1, 2}, {3, 4}});
    int scalar1 = 2;
    mat14 *= scalar1;
    assert(mat14 == Matrix<int>({{2, 4}, {6, 8}}));

    // Test 16.2: Scalar multiplication for a 3x3 matrix
    Matrix<int> mat15({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    int scalar2 = 3;
    mat15 *= scalar2;
    assert(mat15 == Matrix<int>({{3, 6, 9}, {12, 15, 18}, {21, 24, 27}}));

    // Test 17: Matrix transpose
    std::cout << "Test 17: Matrix transpose..." << std::endl;
    Matrix<int> mat16(3, 2);
    mat16(0, 0) = 1;
    mat16(0, 1) = 2;
    mat16(1, 0) = 3;
    mat16(1, 1) = 4;
    mat16(2, 0) = 5;
    mat16(2, 1) = 6;

    Matrix<int> transposed16 = mat16.transpose();

    assert(transposed16.getRowSize() == 2);
    assert(transposed16.getColSize() == 3);
    assert(transposed16(0, 0) == 1);
    assert(transposed16(0, 1) == 3);
    assert(transposed16(0, 2) == 5);
    assert(transposed16(1, 0) == 2);
    assert(transposed16(1, 1) == 4);
    assert(transposed16(1, 2) == 6);

    std::cout << "Test 18: Inequality operator..." << std::endl;

    Matrix<int> mat17({{1, 2}, {3, 4}});
    Matrix<int> mat18({{1, 2}, {3, 4}});
    Matrix<int> mat19({{1, 2, 3}, {4, 5, 6}});

    assert(mat17 != mat19);    // Matrices have different sizes
    assert(!(mat17 != mat18)); // Matrices are equal

    // Test 19: Determinant of a square matrix
    std::cout << "Test 19: Determinant of a square matrix..." << std::endl;
    Matrix<int> mat20({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    assert(mat20.determinant() == 0); // The determinant of this matrix is 0

    // Test 20: Inverse of a square matrix
    std::cout << "Test 20: Inverse of a square matrix..." << std::endl;
    Matrix<double> mat21({{1.0, 2.0, 1.0}, {0.0, 1.0, 2.0}, {3.0, 0.0, 1.0}});
    Matrix<double> inverse21 = mat21.inverse();
    Matrix<double> identity21 = mat21 * inverse21;
    const double tolerance = 1e-6;

    for (int i = 0; i < mat21.getRowSize(); i++)
    {
        for (int j = 0; j < mat21.getColSize(); j++)
        {
            if (i == j)
            {
                assert(std::abs(identity21(i, j) - 1.0) < tolerance);
            }
            else
            {
                assert(std::abs(identity21(i, j) - 0.0) < tolerance);
            }
        }
    }

    // Test 21: Minor matrix
    std::cout << "Test 21: Minor matrix..." << std::endl;
    Matrix<int> mat22({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    Matrix<int> minor22 = mat22.getMinor(0, 1);
    assert(minor22 == Matrix<int>({{4, 6}, {7, 9}}));

    // Test 22: Full rank of a matrix
    std::cout << "Test 22: Full rank of a matrix..." << std::endl;
    Matrix<int> mat23({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    assert(!mat23.isFullRank());

    Matrix<int> mat24({{1, 0, 0}, {0, 2, 0}, {0, 0, 3}});
    assert(mat24.isFullRank());

    // Test 23: Get rank of a matrix
    std::cout << "Test 23: Get rank of a matrix..." << std::endl;
    Matrix<int> mat25({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    assert(mat25.getRank() == 2);

    Matrix<int> mat26({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}});
    assert(mat26.getRank() == 2);

    Matrix<int> mat27({{1, 0, 0}, {0, 2, 0}, {0, 0, 3}});
    assert(mat27.getRank() == 3);

    Matrix<int> mat28({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
    assert(mat28.getRank() == 0);

    // Test 24: Swap rows in a matrix
    std::cout << "Test 24: Swap rows in a matrix..." << std::endl;
    Matrix<int> mat29({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    Matrix<int> expectedMat({{4, 5, 6}, {1, 2, 3}, {7, 8, 9}});
    mat29.swapRows(0, 1);
    assert(mat29 == expectedMat);

    // Test 25: Upper Triangular
    std::cout << "Test 25: Upper Triangular..." << std::endl;
    Matrix<double> matUpperTriangular({{1.0, 2.0, 3.0},
                                       {0.0, 4.0, 5.0},
                                       {0.0, 0.0, 6.0}});
    Matrix<double> upperTriangular = matUpperTriangular.toUpperTriangular();
    assert(upperTriangular(1, 0) == 0.0);
    assert(upperTriangular(2, 0) == 0.0);
    assert(upperTriangular(2, 1) == 0.0);

    // Test 26: Lower Triangular
    std::cout << "Test 26: Lower Triangular..." << std::endl;

    Matrix<double> matLowerTriangular({{1.0, 0.0, 0.0},
                                       {2.0, 3.0, 0.0},
                                       {4.0, 5.0, 6.0}});
    Matrix<double> lowerTriangular = matLowerTriangular.toLowerTriangular();
    assert(lowerTriangular(0, 1) == 0.0);
    assert(lowerTriangular(0, 2) == 0.0);
    assert(lowerTriangular(1, 2) == 0.0);

    // Test 27: Diagonal
    std::cout << "Test 27: Diagonal..." << std::endl;

    Matrix<double> matDiagonal({{1.0, 0.0, 0.0},
                                {0.0, 2.0, 0.0},
                                {0.0, 0.0, 3.0}});
    Matrix<double> diagonal = matDiagonal.toDiagonal();
    assert(diagonal(0, 1) == 0.0);
    assert(diagonal(0, 2) == 0.0);
    assert(diagonal(1, 2) == 0.0);

    // Test 28: Constructing a zeros matrix
    std::cout << "Test 28: Constructing a zeros matrix..." << std::endl;
    Matrix<int> zeros = Matrix<int>::zeros(2, 3);
    assert(zeros(0, 0) == 0);
    assert(zeros(0, 2) == 0);

    // Test 29: Constructing a ones matrix
    std::cout << "Test 29: Constructing a ones matrix..." << std::endl;
    Matrix<int> ones = Matrix<int>::ones(2, 3);
    assert(ones(0, 0) == 1);
    assert(ones(0, 2) == 1);

    // Test 30: Constructing a ones matrix
    std::cout << "Test 30: Constructing an identity matrix..." << std::endl;
    Matrix<int> identity1 = Matrix<int>::identity(3, 3);
    Matrix<int> identity2 = Matrix<int>::identity(5);
    assert(identity1(0, 0) == 1);
    assert(identity1(1, 1) == 1);
    assert(identity1(0, 1) == 0);
    assert(identity2(2, 2) == 1);
    assert(identity2(4, 4) == 1);
    assert(identity2(1, 2) == 0);

    // Test case 31: Test Matrix::random
    std::cout << "Test 31: Constructing a random matrix..." << std::endl;
    Matrix<float> randomMatrix = Matrix<float>::random(0.0f, 10.0f, 10, 10);
    float val1 = randomMatrix(0, 0);
    float val2 = randomMatrix(1, 1);
    assert(val1 >= 0.0f && val1 <= 10.0f);
    assert(val2 >= 0.0f && val2 <= 10.0f);

    // Test case 32: Test Matrix::min
    std::cout << "Test 32: Find minimum value of matrix..." << std::endl;
    Matrix<int> minMatrix(2, 3, {5, 2, 7, 1, 9, 4});
    int minVal = minMatrix.getMin();
    assert(minVal == 1);

    // Test case 33: Test Matrix::min
    std::cout << "Test 33: Find maxmimum value of matrix..." << std::endl;
    Matrix<int> maxMatrix(2, 3, {5, 2, 7, 1, 9, 4});
    int maxVal = maxMatrix.getMax();
    assert(maxVal == 9);

    // Test case 34: Test Matrix::mean
    std::cout << "Test 34: Find mean value of int and float matrix..." << std::endl;
    Matrix<float> floatMatrix(2, 3, {5.0f, -2.0f, -7.0f, 1.0f, 9.0f, 4.0f});
    int meanValInt = minMatrix.mean();
    float meanValFloat = floatMatrix.mean();
    assert(meanValInt == (5 + 2 + 7 + 1 + 9 + 4) / 6);
    assert(meanValFloat == (5.0f - 2.0f - 7.0f + 1.0f + 9.0f + 4.0f) / 6.0f);

    // Test case 34: Test Matrix::abs
    std::cout << "Test 35: Find absolute value of a matrix..." << std::endl;
    Matrix<float> mat30(3, 2, {5.0f, -2.0f, -7.0f, 1.0f, 9.0f, 4.0f});
    Matrix<float> absMatrix = mat30.abs();
    assert(absMatrix(0, 0) == 5.0f);
    assert(absMatrix(1, 1) == 1.0f);
    assert(absMatrix(0, 1) == 2.0f);
    assert(absMatrix(1, 0) == 7.0f);

    // Test case 34: Test Matrix::swapVals
    std::cout << "Test 35: Swap values of a matrix..." << std::endl;
    Matrix<int> swapMatrix(2, 2, {5, 2, 7, 1});
    swapMatrix.swapVals(5, 1);
    assert(swapMatrix(0, 0) == 1);
    assert(swapMatrix(1, 0) == 7);
    assert(swapMatrix(0, 1) == 2);
    assert(swapMatrix(1, 1) == 5);

    // Test case 35: Test Matrix::sort (ascending order)
    std::cout << "Test 35: Sort values of a matrix in ascending order..." << std::endl;
    Matrix<int> sortMatrix(3, 3, {5, 2, 7, 1, 4, 6, 3, 9, 8});
    sortMatrix.sort();
    assert(sortMatrix(0, 0) == 1);
    assert(sortMatrix(0, 1) == 2);
    assert(sortMatrix(0, 2) == 3);
    assert(sortMatrix(1, 0) == 4);
    assert(sortMatrix(1, 1) == 5);
    assert(sortMatrix(1, 2) == 6);
    assert(sortMatrix(2, 0) == 7);
    assert(sortMatrix(2, 1) == 8);
    assert(sortMatrix(2, 2) == 9);

    // Test case 36: Test Matrix::reverseSort (descending order)
    std::cout << "Test 36: Sort values of a matrix in descending order..." << std::endl;
    Matrix<int> reverseSortMatrix(3, 3, {5, 2, 7, 1, 4, 6, 3, 9, 8});
    reverseSortMatrix.reverseSort();
    assert(reverseSortMatrix(0, 0) == 9);
    assert(reverseSortMatrix(0, 1) == 8);
    assert(reverseSortMatrix(0, 2) == 7);
    assert(reverseSortMatrix(1, 0) == 6);
    assert(reverseSortMatrix(1, 1) == 5);
    assert(reverseSortMatrix(1, 2) == 4);
    assert(reverseSortMatrix(2, 0) == 3);
    assert(reverseSortMatrix(2, 1) == 2);
    assert(reverseSortMatrix(2, 2) == 1);

    // Test case 36: Test Matrix::median for even-sized matrix
    std::cout << "Test 36: Median of a matrix (even-sized)..." << std::endl;
    Matrix<int> evenSizedMatrix(2, 2, {3, 2, 1, 4});
    int medianEven = evenSizedMatrix.median();
    Matrix<float> evenSizedMatrix2(2, 2, {3, 2, 1, 4});
    float medianEven2 = evenSizedMatrix2.median();
    assert(medianEven == 2);
    assert(medianEven2 == 2.5);

    // Test case 37: Test Matrix::median for odd-sized matrix
    std::cout << "Test 37: Median of a matrix (odd-sized)..." << std::endl;
    Matrix<int> oddSizedMatrix(3, 3, {7, 1, 5, 2, 4, 6, 3, 9, 8});
    int medianOdd = oddSizedMatrix.median();
    assert(medianOdd == 5);

    // Test case 38: Test Matrix::hstack
    std::cout << "Test 38: Horizontal Stack of Matrices..." << std::endl;
    Matrix<int> matrixA(2, 2, {1, 2, 3, 4});
    Matrix<int> matrixB(2, 3, {5, 6, 7, 8, 9, 10});
    Matrix<int> hstackResult = matrixA.hstack(matrixB);
    assert(hstackResult.getRowSize() == 2);
    assert(hstackResult.getColSize() == 5);
    assert(hstackResult(0, 2) == 5);
    assert(hstackResult(1, 4) == 10);

    // Test case 39: Test Matrix::vstack
    std::cout << "Test 39: Vertical Stack of Matrices..." << std::endl;
    Matrix<int> matrixC(2, 2, {1, 2, 3, 4});
    Matrix<int> matrixD(3, 2, {5, 6, 7, 8, 9, 10});
    Matrix<int> vstackResult = matrixC.vstack(matrixD);
    assert(vstackResult.getRowSize() == 5);
    assert(vstackResult.getColSize() == 2);
    assert(vstackResult(3, 1) == 8);

    // Test case 40: Test Matrix::reshape
    std::cout << "Test 40: Reshaping matrices..." << std::endl;
    Matrix<int> A(2, 3, {1, 2, 3, 4, 5, 6});
    A.reshape(3, 2);
    assert(A.getRowSize() == 3);
    assert(A.getColSize() == 2);
    assert(A(0, 0) == 1);
    assert(A(0, 1) == 2);
    assert(A(1, 0) == 3);
    assert(A(1, 1) == 4);
    assert(A(2, 0) == 5);
    assert(A(2, 1) == 6);

    Matrix<float> f(3, 6, {1.1f, 2.2f, 3.3f, 4.4f, 5.5f, 6.6f, 7.7f, 8.8f, 9.9f, 10.1f, 11.1f, 12.2f, 4.4f, 5.5f, 6.6f, 7.7f, 8.8f, 9.9f});
    f.reshape(6, 3);
    assert(f.getRowSize() == 6);
    assert(f.getColSize() == 3);
    assert(f(0, 0) == 1.1f);
    assert(f(0, 1) == 2.2f);
    assert(f(0, 2) == 3.3f);
    assert(f(1, 0) == 4.4f);
    assert(f(1, 1) == 5.5f);
    assert(f(1, 2) == 6.6f);
    assert(f(2, 0) == 7.7f);
    assert(f(2, 1) == 8.8f);
    assert(f(2, 2) == 9.9f);
    assert(f(3, 0) == 10.1f);
    assert(f(3, 1) == 11.1f);
    assert(f(3, 2) == 12.2f);

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
