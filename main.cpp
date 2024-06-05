#include <iostream>
#include "Matrix.h"
#include <assert.h>
#include "LA.h"
#include "ComplexNum.h"

int main() {
    // Complex number test asserts
    ComplexNum z(3, 2);
    assert(z.getConjugate().getImagPart() == -2 );
    assert ((2 + z).getRealPart() == 5);
    assert((z * -7).getImagPart() == -14);
    assert((2 * z).getRealPart() == 6);
    assert((z * 12).getRealPart() == 36);
    assert((z * 0).getRealPart() == 0);
    assert(z.getMagnitude() == sqrt(13));
    ComplexNum w(2, 3);
    assert((z*w).getImagPart() == 13);
    assert((z*w).getRealPart() == 0);

    // Matrix test asserts
    Matrix A(2, 2);
    A(0,1) = ComplexNum(1, 0);
    assert(A(0, 1) == ComplexNum(1, 0));
    Matrix B(2, 2);
    B(0, 1) = ComplexNum(3, 4);
    Matrix C = A + B;
    assert(C(0, 1) == ComplexNum(4, 4));

    Matrix CTranspose = conjTranspose(&C);
    assert(CTranspose(1, 0) == ComplexNum(4, -4));
    assert(CTranspose(0, 1) == ComplexNum(0, 0));

    Matrix matrixToNorm(3, 2);
    matrixToNorm(0, 0) = 3;
    matrixToNorm(1 ,0) = 4;
    matrixToNorm(2, 0) = 2;
    matrixToNorm(0, 1) = 1;
    matrixToNorm(1, 1) = ComplexNum(0, 2);
    matrixToNorm(2, 1) = 1;
    std::cout << "Frobenius norm: " << frobeniusNorm(&matrixToNorm) << std::endl;



    // Gram-Schmidt test asserts
    Matrix D(3,3);
    D(0, 0) = 1;
    D(0,1) = 8;

    D(1,0) = 2;
    D(1,1) = 1;
    D(2, 1) = -6;
    D(2,2) = 1;
    auto meme = D[1];
    auto G = GramSchmidt(D);
    std::cout << ComplexNum(-4, -4) << std::endl;
    std::cout << D << std::endl;
}
