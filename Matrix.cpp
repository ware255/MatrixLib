#include <cstdlib>
#include <cassert>
#include "Matrix.hpp"

template<typename T>
Matrix<T>::Matrix() {
    row = 1;
    col = 1;
    data = (T*)malloc(sizeof(T) * row * col);
    assert(data != NULL);
}

template<typename T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols) {
    assert(rows > 0 && cols > 0);
    row = rows;
    col = cols;
    data = (T*)malloc(sizeof(T) * row * col);
    assert(data != NULL);
}

template<typename T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const T init_num) {
    assert(rows > 0 && cols > 0);
    row = rows;
    col = cols;
    data = (T*)malloc(sizeof(T) * row * col);
    assert(data != NULL);
    for (unsigned int i = 0; i < row * col; ++i)
        data[i] = init_num;
}

template<typename T>
Matrix<T>::~Matrix() {
    free(data);
    data = NULL;
}

template<typename T>
T& Matrix<T>::operator()(unsigned int r, unsigned int c) {
    assert(r < row && c < col);
    return data[c * row + r]; //array[n * M + m]; == array(m,n)
}

template<typename T>
const T& Matrix<T>::operator()(unsigned int r, unsigned int c) const {
    assert(r < row && c < col);
    return data[c * row + r];
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
    if (this != &other) {
        this->row = other.row;
        this->col = other.col;
        free(this->data);
        this->data = (T*)malloc(sizeof(T) * this->row * this->col);
        assert(this->data != NULL);
        memcpy(this->data, other.data, sizeof(T) * this->row * this->col);
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {
    assert(this->row == other.row && this->col == other.col);
    unsigned int i, j;
    for (i = 0; i < this->row; ++i)
        for (j = 0; j < this->col; ++j)
            this->data[j * this->row + i] += other.data[j * other.row + i];
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other) {
    assert(this->row == other.row && this->col == other.col);
    unsigned int i, j;
    for (i = 0; i < this->row; ++i)
        for (j = 0; j < this->col; ++j)
            this->data[j * this->row + i] -= other.data[j * other.row + i];
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T other) {
    unsigned int i, j;
    for (i = 0; i < this->row; ++i)
        for (j = 0; j < this->col; ++j)
            this->data[j * this->row + i] *= other;
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other) {
    assert(this->col == other.row);

    Matrix<T> res(this->row, other.col);
    unsigned int i, j, k;

    for (i = 0; i < this->row; ++i)
        for (j = 0; j < other.col; ++j) {
            res.data[j * res.row + i] = 0;
            for (k = 0; k < this->col; k++)
                res.data[j * res.row + i] += this->data[k * this->row + i] * other.data[j * other.row + k];
        }
    *this = res;

    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    Matrix<T> res(this->row, this->col);
    res = *this;
    res += other;
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    Matrix<T> res(this->row, this->col);
    res = *this;
    res -= other;
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T other) const {
    Matrix<T> res(this->row, this->col);
    res = *this;
    res *= other;
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
    Matrix<T> res(this->row, this->col);
    res = *this;
    res *= other;
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator%(const T mod) {
    Matrix<T> res(this->row, this->col);
    for (unsigned int i = 0; i < this->row; ++i)
        for (unsigned int j = 0; j < this->col; ++j)
            res.data[j * this->row + i] = this->data[i * this->row + j] % mod;
    return res;
}


template<typename T>
T Matrix<T>::determinant() {
    assert(this->row == this->col);

    double *sweep = (double*)malloc(sizeof(double) * this->row * this->col);
    assert(sweep != NULL);

    double det, max, tmp, a;
    double t_det, multiple, sweep_ik;

    unsigned int i, j, k, N, max_i;

    N = this->row;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            sweep[j * this->row + i] = (double)this->data[j * this->row + i];

    multiple = 1.0;

    for (k = 0; k < N; k++) {
        max = (sweep[k * this->row + k] < 0) ?
            -sweep[k * this->row + k] : sweep[k * this->row + k];
        max_i = k;

        for (i = k + 1; i < N; i++) {
            sweep_ik = (sweep[k * this->row + i] < 0) ?
                -sweep[k * this->row + i] : sweep[k * this->row + i];
            if (sweep_ik > max) {
                max = sweep_ik;
                max_i = i;
            }
        }

        if (((sweep[k * this->row + max_i] < 0) ?
                -sweep[k * this->row + max_i] : sweep[k * this->row + max_i]) <= 1e-9) {
            multiple = 0.0;
            break;
        }

        if (k != max_i) {
            for (j = 0; j < N; j++) {
                tmp = sweep[j * this->row + max_i];
                sweep[j * this->row + max_i] = sweep[j * this->row + k];
                sweep[j * this->row + k] = tmp;
            }

            multiple *= -1.0;
        }

        a = 1.0 / sweep[k * this->row + k];

        for (j = 0; j < N; j++)
            sweep[j * this->row + k] *= a;

        multiple *= a;

        for (i = 0; i < N; i++) {
            if (i == k) continue;

            a = -sweep[k * this->row + i];

            for (j = 0; j < N; j++)
                sweep[j * this->row + i] += sweep[j * this->row + k] * a; 
        }
    }

    if (multiple == 0) {
        det = 0;
    } else {
        t_det = 1;
        det = t_det / multiple;
    }

    free(sweep);

    return static_cast<T>(det);
}

template<typename T>
Matrix<T> Matrix<T>::transpose() {
    Matrix<T> res(this->col, this->row);
    unsigned int i, j;
    for (i = 0; i < this->row; ++i)
        for (j = 0; j < this->col; ++j)
            res.data[i * this->col + j] = this->data[j * this->row + i];
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::inverse() {
    assert(this->row == this->col);

    Matrix<T> A(this->row, this->col);
    Matrix<T> res(this->row, this->col);
    unsigned int i, j, k, N;

    A = *this;
    N = this->row;
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j)
            res.data[j * this->row + i] = (i == j);

    for (i = 0; i < N; ++i) {
        T buf = 1 / A.data[i * this->row + i];
        for (j = 0; j < N; ++j) {
            A.data[j * this->row + i] *= buf;
            res.data[j * this->row + i] *= buf;
        }
        for (j = 0; j < N; ++j)
            if (i != j) {
                buf = A.data[i * this->row + j];
                for (k = 0; k < N; ++k) {
                    A.data[k * this->row + j] -= A.data[k * this->row + i] * buf;
                    res.data[k * this->row + j] -= res.data[k * this->row + i] * buf;
                }
            }
    }

    return res;
}

template<typename T>
bool Matrix<T>::is_identity_matrix() {
    bool flag = true;
    unsigned int i, j;
    for (i = 0; i < this->row; ++i)
        for (j = 0; j < this->col; ++j)
            if (i == j)
                if (this->data[j * this->row + i] == 0)
                    flag = false;
            else
                if (this->data[j * this->row + i] == 1)
                    flag = false;
    return flag;
}

template<typename T>
bool Matrix<T>::is_zero_matrix() {
    unsigned int i, j;
    for (i = 0; i < this->row; ++i)
        for (j = 0; j < this->col; ++j)
            if (this->data[j * this->row + i] != 0)
                return false;
    return true;
}
