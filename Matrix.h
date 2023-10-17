#ifndef MATRIX_IS_INCLUDED
#define MATRIX_IS_INCLUDED

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <vector>

template <class T>
class Matrix
{
private:
    int rows, cols;
    T *pMatrix;

public:
    Matrix()
    {
        this->rows = 0;
        this->cols = 0;
        this->pMatrix = new T[rows * cols];
    }

    Matrix(int rows, int cols)
    {
        this->rows = rows;
        this->cols = cols;
        this->pMatrix = new T[rows * cols];
    }

    Matrix(const Matrix &m)
    {
        this->rows = m.rows;
        this->cols = m.cols;
        this->pMatrix = new T[rows * cols];

        for (int i = 0; i < rows * cols; i++)
        {
            this->pMatrix[i] = m.pMatrix[i];
        }
    }

    Matrix(const std::vector<std::vector<T>> &values)
    {
        this->rows = values.size();
        if (this->rows == 0)
        {
            this->cols = 0;
            this->pMatrix = nullptr;
            return;
        }
        this->cols = values[0].size();

        // Allocate memory and copy values
        this->pMatrix = new T[rows * cols];
        for (int i = 0; i < rows; i++)
        {
            if (values[i].size() != cols)
            {
                throw std::invalid_argument("All rows must have the same number of columns.");
            }
            for (int j = 0; j < cols; j++)
            {
                this->pMatrix[i * cols + j] = values[i][j];
            }
        }
    }

    Matrix(int rowSize, int colSize, const std::vector<T> &values)
    {
        if (rowSize * colSize != values.size())
        {
            throw std::invalid_argument("Size of the values vector must match the specified matrix dimensions.");
        }

        this->rows = rowSize; // Update rows
        this->cols = colSize; // Update cols

        // Allocate memory and copy values
        this->pMatrix = new T[rows * cols];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                this->pMatrix[i * cols + j] = values[i * colSize + j];
            }
        }
    }

    ~Matrix()
    {
        if (nullptr != pMatrix)
        {
            delete[] pMatrix;
        }
    }

    T &operator[](int index)
    {
        if (index < 0 || index >= rows * cols)
        {
            throw std::out_of_range("Index out of range.");
        }
        return pMatrix[index];
    }

    T &operator()(int row, int col)
    {
        if (row < 0 || row >= rows || col < 0 || col >= cols)
        {
            throw std::out_of_range("Index out of range.");
        }
        return pMatrix[row * cols + col];
    }

    // Const version for reading elements from a const object
    const T &operator()(int row, int col) const
    {
        if (row < 0 || row >= rows || col < 0 || col >= cols)
        {
            throw std::out_of_range("Index out of range.");
        }
        return pMatrix[row * cols + col];
    }

    /*
        Getters and setters
    */
    int getRowSize() const
    {
        return rows;
    }

    int getColSize() const
    {
        return cols;
    }

    int getSize() const
    {
        return rows * cols;
    }

    void setRowSize(const T &value)
    {
        rows = value;

        if (nullptr != pMatrix)
        {
            delete[] pMatrix;
        }

        pMatrix = new T[rows * cols];
    }

    void setColSize(const T &value)
    {
        cols = value;

        if (nullptr != pMatrix)
        {
            delete[] pMatrix;
        }

        pMatrix = new T[rows * cols];
    }

    /*
    Arithematic operators
    */

    Matrix<T> &operator=(const Matrix<T> &m)
    {
        if (this == &m)
            return *this;

        if (rows != m.rows || cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be the same for assignment.");
        }

        for (int i = 0; i < rows * cols; i++)
        {
            pMatrix[i] = m.pMatrix[i];
        }

        return *this;
    }

    bool operator==(const Matrix<T> &m) const
    {
        if (rows != m.rows || cols != m.cols)
            return false;

        for (int i = 0; i < rows * cols; i++)
        {
            if (pMatrix[i] != m.pMatrix[i])
                return false;
        }

        return true;
    }

    bool operator!=(const Matrix<T> &m) const
    {
        return !(*this == m);
    }

    Matrix<T> operator+(const Matrix<T> m) const
    {
        if (rows != m.rows || cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be equal for matrix addition.");
        }

        Matrix<T> result(rows, cols);

        for (int i = 0; i < rows * cols; i++)
        {
            result.pMatrix[i] = pMatrix[i] + m.pMatrix[i];
        }

        return result;
    }

    Matrix<T> operator-(const Matrix<T> m) const
    {
        if (rows != m.rows || cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be equal for matrix subtraction.");
        }

        Matrix<T> result(rows, cols);

        for (int i = 0; i < rows * cols; i++)
        {
            result.pMatrix[i] = pMatrix[i] - m.pMatrix[i];
        }

        return result;
    }

    Matrix<T> operator*(const Matrix<T> &m) const
    {
        if (cols != m.rows)
        {
            throw std::invalid_argument("Number of columns in the first matrix must be equal to the number of rows in the second matrix for multiplication.");
        }

        Matrix<T> result(rows, m.cols);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < m.cols; j++)
            {
                result.pMatrix[i * result.cols + j] = 0;
                for (int k = 0; k < cols; k++)
                {
                    result.pMatrix[i * result.cols + j] += pMatrix[i * cols + k] * m.pMatrix[k * m.cols + j];
                }
            }
        }

        return result;
    }

    Matrix<T> operator*(const T &scalar) const
    {
        Matrix<T> result(rows, cols);

        for (int i = 0; i < rows * cols; i++)
        {
            result.pMatrix[i] = pMatrix[i] * scalar;
        }

        return result;
    }

    Matrix<T> operator/(const T &scalar) const
    {
        Matrix<T> result(rows, cols);

        for (int i = 0; i < rows * cols; i++)
        {
            result.pMatrix[i] = pMatrix[i] / scalar;
        }

        return result;
    }

    Matrix<T> &operator+=(const Matrix<T> &m)
    {
        if (rows != m.rows || cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be the same for addition assignment.");
        }

        for (int i = 0; i < rows * cols; i++)
        {
            pMatrix[i] += m.pMatrix[i];
        }

        return *this;
    }

    Matrix<T> &operator-=(const Matrix<T> &m)
    {
        if (rows != m.rows || cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be the same for subtraction assignment.");
        }

        for (int i = 0; i < rows * cols; i++)
        {
            pMatrix[i] -= m.pMatrix[i];
        }

        return *this;
    }

    Matrix<T> &operator*=(const Matrix<T> &m)
    {
        if (this->getColSize() != m.getRowSize())
        {
            throw std::invalid_argument("Number of columns in the first matrix must be equal to the number of rows in the second matrix for multiplication assignment.");
        }

        Matrix<T> result(this->getRowSize(), m.getColSize());

        for (int i = 0; i < this->getRowSize(); i++)
        {
            for (int j = 0; j < m.getColSize(); j++)
            {
                result(i, j) = 0;
                for (int k = 0; k < this->getColSize(); k++)
                {
                    result(i, j) += (*this)(i, k) * m(k, j);
                }
            }
        }

        this->rows = result.getRowSize();
        this->cols = result.getColSize();
        delete[] this->pMatrix;
        this->pMatrix = new T[this->rows * this->cols];
        for (int i = 0; i < this->rows; i++)
        {
            for (int j = 0; j < this->cols; j++)
            {
                (*this)(i, j) = result(i, j);
            }
        }

        return *this;
    }

    Matrix<T> &operator*=(const T &scalar)
    {
        for (int i = 0; i < rows * cols; i++)
        {
            pMatrix[i] *= scalar;
        }
        return *this;
    }

    /*
        Linear algebra functions
    */

    Matrix<T> transpose() const
    {
        Matrix<T> result(cols, rows);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                result(j, i) = (*this)(i, j); // Swap rows and columns
            }
        }

        return result;
    }

    T determinant() const
    {
        if (this->getRowSize() != this->getColSize())
        {
            throw std::invalid_argument("Matrix must be square for determinant calculation.");
        }

        int size = this->getRowSize();
        if (size == 1)
        {
            return (*this)(0, 0);
        }
        if (size == 2)
        {
            return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
        }

        T det = 0;
        for (int j = 0; j < size; j++)
        {
            Matrix<T> minor = this->getMinor(0, j);
            det += ((*this)(0, j) * minor.determinant() * ((j % 2 == 0) ? 1 : -1));
        }
        return det;
    }

    Matrix<T> inverse() const
    {
        if (this->getRowSize() != this->getColSize())
        {
            throw std::invalid_argument("Matrix must be square for matrix inversion.");
        }

        T det = this->determinant();
        if (det == 0)
        {
            throw std::invalid_argument("Matrix is singular; it doesn't have an inverse.");
        }

        int size = this->getRowSize();
        Matrix<T> adjugate(size, size);

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                Matrix<T> minor = this->getMinor(i, j);
                T cofactor = minor.determinant() * (((i + j) % 2 == 0) ? 1 : -1);
                adjugate(j, i) = cofactor; // Swap indices here
            }
        }

        Matrix<T> inverseMatrix = adjugate / det;
        return inverseMatrix;
    }

    Matrix<T> getMinor(int rowToRemove, int colToRemove) const
    {
        if (rowToRemove < 0 || rowToRemove >= this->getRowSize() || colToRemove < 0 || colToRemove >= this->getColSize())
        {
            throw std::out_of_range("Row or column index out of range.");
        }

        int minorRows = this->getRowSize() - 1;
        int minorCols = this->getColSize() - 1;

        Matrix<T> minor(minorRows, minorCols);

        for (int i = 0, minorRow = 0; i < this->getRowSize(); i++)
        {
            if (i == rowToRemove)
            {
                continue;
            }

            for (int j = 0, minorCol = 0; j < this->getColSize(); j++)
            {
                if (j == colToRemove)
                {
                    continue;
                }

                minor(minorRow, minorCol) = (*this)(i, j);
                minorCol++;
            }

            minorRow++;
        }

        return minor;
    }

    void swapRows(int row1, int row2)
    {
        if (row1 == row2 || row1 < 0 || row1 >= rows || row2 < 0 || row2 >= rows)
        {
            throw std::out_of_range("Row indices are out of range or identical.");
        }

        for (int col = 0; col < cols; col++)
        {
            std::swap((*this)(row1, col), (*this)(row2, col));
        }
    }

    int getRank() const
    {
        Matrix<T> copy(*this);
        int rank = 0;
        int numRows = copy.getRowSize();
        int numCols = copy.getColSize();

        for (int row = 0; row < numRows; row++)
        {
            int leadingColumn = -1;
            for (int col = 0; col < numCols; col++)
            {
                if (copy(row, col) != 0)
                {
                    leadingColumn = col;
                    break;
                }
            }

            if (leadingColumn != -1)
            {
                rank++;

                T leadingElement = copy(row, leadingColumn);
                for (int col = leadingColumn; col < numCols; col++)
                {
                    copy(row, col) /= leadingElement;
                }

                for (int r = 0; r < numRows; r++)
                {
                    if (r != row && copy(r, leadingColumn) != 0)
                    {
                        T factor = -copy(r, leadingColumn);
                        for (int col = leadingColumn; col < numCols; col++)
                        {
                            copy(r, col) += factor * copy(row, col);
                        }
                    }
                }
            }
        }

        return rank;
    }

    bool isFullRank() const
    {
        int rank = getRank();
        int minDim = std::min(rows, cols);

        return rank == minDim;
    }

    Matrix<T> toUpperTriangular() const
    {
        if (this->getRowSize() != this->getColSize())
        {
            throw std::invalid_argument("Matrix must be square for triangular conversion.");
        }

        Matrix<T> result(*this);

        for (int i = 0; i < this->getRowSize(); i++)
        {
            for (int j = 0; j < i; j++)
            {
                result(i, j) = 0;
            }
        }

        return result;
    }

    Matrix<T> toLowerTriangular() const
    {
        if (this->getRowSize() != this->getColSize())
        {
            throw std::invalid_argument("Matrix must be square for triangular conversion.");
        }

        Matrix<T> result(*this);

        for (int i = 0; i < this->getRowSize(); i++)
        {
            for (int j = i + 1; j < this->getColSize(); j++)
            {
                result(i, j) = 0;
            }
        }

        return result;
    }

    Matrix<T> toDiagonal() const
    {
        if (this->getRowSize() != this->getColSize())
        {
            throw std::invalid_argument("Matrix must be square for diagonal conversion.");
        }

        Matrix<T> result(*this);

        for (int i = 0; i < this->getRowSize(); i++)
        {
            for (int j = 0; j < this->getColSize(); j++)
            {
                if (i != j)
                {
                    result(i, j) = 0;
                }
            }
        }

        return result;
    }

    // void LUDecomp(Matrix<T> &L, Matrix<T> &U)
    // {
    //     if (this->getRowSize() != this->getColSize())
    //     {
    //         throw std::invalid_argument("Matrix must be square for LU decomposition.");
    //     }

    //     int size = this->getRowSize();

    //     L = Matrix<T>(size, size);
    //     U = Matrix<T>(size, size);

    //     for (int i = 0; i < size; i++)
    //     {
    //         // U matrix
    //         for (int j = i; j < size; j++)
    //         {
    //             T sum = 0;
    //             for (int k = 0; k < i; k++)
    //             {
    //                 sum += L(i, k) * U(k, j);
    //             }
    //             U(i, j) = (*this)(i, j) - sum;
    //         }

    //         // L matrix
    //         for (int j = i; j < size; j++)
    //         {
    //             if (i == j)
    //                 L(i, j) = 1;
    //             else
    //             {
    //                 T sum = 0;
    //                 for (int k = 0; k < i; k++)
    //                 {
    //                     sum += L(j, k) * U(k, i);
    //                 }
    //                 L(j, i) = ((*this)(j, i) - sum) / U(i, i);
    //             }
    //         }
    //     }
    // }

    // Matrix<T> choleskyDecomp() const
    // {
    //     if (*this != this->transpose())
    //     {
    //         throw std::invalid_argument("Matrix must be symmetric for Cholesky decomposition.");
    //     }

    //     Matrix result = *this;

    //     for (int i = 0; i < this->rows; i++)
    //     {
    //         for (int j = 0; j < this->cols; j++)
    //         {
    //             if (i == j)
    //             {
    //                 for (int k = 0; k < j; k++)
    //                 {
    //                     result(i, j) -= pow(result(i, k), 2);
    //                 }

    //                 result(i, j) = sqrt(result(i, j));
    //             }

    //             else if (i > j)
    //             {
    //                 for (int k = 0; k < j; k++)
    //                 {
    //                     result(i, j) -= result(i, k) * result(j, k);
    //                 }

    //                 result(i, j) /= result(j, j);
    //             }

    //             else
    //             {
    //                 result(i, j) = T();
    //             }
    //         }
    //     }

    //     return result;
    // }

    /*
        Stream operators
    */

    template <class U>
    friend std::ostream &operator<<(std::ostream &os, const Matrix<U> &matrix);

    // Friend stream input operator
    template <class U>
    friend std::istream &operator>>(std::istream &is, Matrix<U> &matrix);

    template <class U>
    friend std::ofstream &operator<<(std::ofstream &ofs, const Matrix<U> &matrix);

    template <class U>
    friend std::ifstream &operator>>(std::ifstream &ifs, Matrix<U> &matrix);
};

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix)
{
    for (int i = 0; i < matrix.rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            os << matrix.pMatrix[i * matrix.cols + j] << ' ';
        }
        os << '\n';
    }
    return os;
}

// Friend stream input operator for Matrix
template <class T>
std::istream &operator>>(std::istream &is, Matrix<T> &matrix)
{
    for (int i = 0; i < matrix.rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            is >> matrix.pMatrix[i * matrix.cols + j];
        }
    }
    return is;
}

template <class T>
std::ofstream &operator<<(std::ofstream &ofs, const Matrix<T> &matrix)
{
    for (int i = 0; i < matrix.rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            ofs << matrix.pMatrix[i * matrix.cols + j] << ' ';
        }
        ofs << '\n';
    }
    return ofs;
}

template <class T>
std::ifstream &operator>>(std::ifstream &ifs, Matrix<T> &matrix)
{
    for (int i = 0; i < matrix.rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            ifs >> matrix.pMatrix[i * matrix.cols + j];
        }
    }
    return ifs;
}

template <class T>
class SquareMatrix : public Matrix<T>
{
private:
    int rows, cols;
    T *pMatrix;

public:
    SquareMatrix(int size) : Matrix<T>(size, size) {}
};

#endif