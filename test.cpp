#include <cstdio>
#include <cmath>
#include "Matrix.hpp"

int main(void) {
    Matrix<int> mat(4,4);
    mat(0,0)=3; mat(0,1)=1; mat(0,2)=1; mat(0,3)=2;
    mat(1,0)=5; mat(1,1)=1; mat(1,2)=3; mat(1,3)=4;
    mat(2,0)=2; mat(2,1)=0; mat(2,2)=1; mat(2,3)=0;
    mat(3,0)=1; mat(3,1)=3; mat(3,2)=2; mat(3,3)=1;

    puts("mat");
    for (unsigned int i = 0; i < mat.row_size(); i++) {
        for (unsigned int j = 0; j < mat.col_size(); j++)
            printf("%5d ", mat(i,j));
        puts("");
    }

    mat = mat * 2;
    puts("mat * 2");
    for (unsigned int i = 0; i < mat.row_size(); i++) {
        for (unsigned int j = 0; j < mat.col_size(); j++)
            printf("%5d ", mat(i,j));
        puts("");
    }

    return 0;
}
