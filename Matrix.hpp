#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

template<typename T>
class Matrix {
    unsigned int row; // 行
    unsigned int col; // 列
    T *data;
public:
    Matrix();
    Matrix(unsigned int rows, unsigned int cols);
    Matrix(unsigned int rows, unsigned int cols, const T init_num);
    ~Matrix();

    unsigned int row_size() const { return row; }
    unsigned int col_size() const { return col; }

    T& operator()(unsigned int r, unsigned int c);
    const T& operator()(unsigned int r, unsigned int c) const;

    Matrix<T>& operator=(const Matrix<T>& other);

    Matrix<T>& operator+=(const Matrix<T>& other);
    Matrix<T>& operator-=(const Matrix<T>& other);
    Matrix<T>& operator*=(const T other);
    Matrix<T>& operator*=(const Matrix<T>& other);
    Matrix<T> operator+(const Matrix<T>& other) const;
    Matrix<T> operator-(const Matrix<T>& other) const;
    Matrix<T> operator*(const T other) const;
    Matrix<T> operator*(const Matrix<T>& other) const;
    Matrix<T> operator%(const T mod);

    T determinant();
    Matrix<T> transpose();
    Matrix<T> inverse();

    bool is_identity_matrix();
    bool is_zero_matrix();
};

#include "Matrix.cpp"

#endif // _MATRIX_HPP_
