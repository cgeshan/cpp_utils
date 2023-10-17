#ifndef MATRIX_IS_INCLUDED
#define MATRIX_IS_INCLUDED

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>

template <class T>
class Matrix
{
private:
    int rows, cols;
    T *pMatrix;

public:
    Matrix(int rows, int cols)
    {
        this->rows = rows;
        this->cols = cols;
        pMatrix = new T[rows * cols];
    }

    Matrix(const Matrix &m)
    {
        this->rows = m.rows;
        this->cols = m.cols;
        pMatrix = new T[rows * cols];

        for (int i = 0; i < rows * cols; i++)
        {
            this->pMatrix[i] = m.pMatrix[i];
        }
    }

    ~Matrix()
    {
        if (nullptr != pMatrix)
        {
            delete[] pMatrix;
        }
    }

    /*
    Arithematic operators
    */

    Matrix<T> &operator=(const Matrix<T> &m)
    {
        if (this == &m)
            return *this;

        if (this->rows != m.rows || this->cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be the same for assignment.");
        }

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            this->pMatrix[i] = m.pMatrix[i];
        }

        return *this;
    }

    bool operator==(const Matrix<T> &m) const
    {
        if (this->rows != m.rows || this->cols != m.cols)
            return false;

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            if (this->pMatrix[i] != m.pMatrix[i])
                return false;
        }

        return true;
    }

    Matrix<T> operator+(const Matrix<T> m) const
    {
        if (this->rows != m.rows || this->cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be equal for matrix addition.");
        }

        Matrix<T> result(this->rows, this->cols);

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            result.pMatrix[i] = this->pMatrix[i] + m.pMatrix[i];
        }

        return result;
    }

    Matrix<T> operator-(const Matrix<T> m) const
    {
        if (this->rows != m.rows || this->cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be equal for matrix subtraction.");
        }

        Matrix<T> result(this->rows, this->cols);

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            result.pMatrix[i] = this->pMatrix[i] - m.pMatrix[i];
        }

        return result;
    }

    Matrix<T> operator*(const Matrix<T> &m) const
    {
        if (this->cols != m.rows)
        {
            throw std::invalid_argument("Number of columns in the first matrix must be equal to the number of rows in the second matrix for multiplication.");
        }

        Matrix<T> result(this->rows, m.cols);

        for (int i = 0; i < this->rows; i++)
        {
            for (int j = 0; j < m.cols; j++)
            {
                result.pMatrix[i * result.cols + j] = 0;
                for (int k = 0; k < this->cols; k++)
                {
                    result.pMatrix[i * result.cols + j] += this->pMatrix[i * this->cols + k] * m.pMatrix[k * m.cols + j];
                }
            }
        }

        return result;
    }

    Matrix<T> operator*(const T &scalar) const
    {
        Matrix<T> result(this->rows, this->cols);

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            result.pMatrix[i] = this->pMatrix[i] * scalar;
        }

        return result;
    }

    Matrix<T> operator/(const T &scalar) const
    {
        Matrix<T> result(this->rows, this->cols);

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            result.pMatrix[i] = this->pMatrix[i] / scalar;
        }

        return result;
    }

    Matrix<T> &operator+=(const Matrix<T> &m)
    {
        if (this->rows != m.rows || this->cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be the same for addition assignment.");
        }

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            this->pMatrix[i] += m.pMatrix[i];
        }

        return *this;
    }

    Matrix<T> &operator-=(const Matrix<T> &m)
    {
        if (this->rows != m.rows || this->cols != m.cols)
        {
            throw std::invalid_argument("Matrix dimensions must be the same for subtraction assignment.");
        }

        for (int i = 0; i < this->rows * this->cols; i++)
        {
            this->pMatrix[i] -= m.pMatrix[i];
        }

        return *this;
    }

    Matrix<T> &operator*=(const Matrix<T> &m)
    {
        if (this->cols != m.rows)
        {
            throw std::invalid_argument("Number of columns in the first matrix must be equal to the number of rows in the second matrix for multiplication assignment.");
        }

        Matrix<T> result(this->rows, m.cols);

        for (int i = 0; i < this->rows; i++)
        {
            for (int j = 0; j < m.cols; j++)
            {
                result.pMatrix[i * result.cols + j] = 0;
                for (int k = 0; k < this->cols; k++)
                {
                    result.pMatrix[i * result.cols + j] += this->pMatrix[i * this->cols + k] * m.pMatrix[k * m.cols + j];
                }
            }
        }

        *this = result;
        return *this;
    }

    Matrix<T> transpose() const
    {
        Matrix<T> result(this->cols, this->rows);

        for (int i = 0; i < this->rows; i++)
        {
            for (int j = 0; j < this->cols; j++)
            {
                result.pMatrix[j * this->rows + i] = this->pMatrix[i * this->cols + j];
            }
        }

        return result;
    }

    T &operator[](int index)
    {
        if (index < 0 || index >= this->rows * this->cols)
        {
            throw std::out_of_range("Index out of range.");
        }
        return this->pMatrix[index];
    }

    /*
        Getters and setters
    */
    int getRowSize()
    {
        return this->rows;
    }

    int getColSize()
    {
        return this->cols;
    }

    int getSize()
    {
        return this->rows * this->cols;
    }

    void setRowSize(const T &value)
    {
        this->rows = value;

        if (nullptr != pMatrix)
        {
            delete[] pMatrix;
        }

        pMatrix = new T[rows * cols];
    }

    void setColSize(const T &value)
    {
        this->cols = value;

        if (nullptr != pMatrix)
        {
            delete[] pMatrix;
        }

        pMatrix = new T[rows * cols];
    }

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