#include <iostream>
#include "Matrix.h"
#include <assert.h>
#include "LA.h"
#include "ComplexNum.h"
#include "Regression.h"

int main() {
    // Complex number test asserts
    ComplexNum z(3, 2);
    assert(z.getConjugate().getImagPart() == -2 );
    assert ((2+z).getRealPart() == 5);
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
    //std::cout << "Frobenius norm: " << frobeniusNorm(&matrixToNorm) << std::endl;



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
    //std::cout << ComplexNum(-4, -4) << std::endl;

    Matrix matrixToInvert(3, 3);
    matrixToInvert(0,0) = ComplexNum(2, 1);
    matrixToInvert(1,0) = ComplexNum(4, 0);
    matrixToInvert(2,0) = ComplexNum(1, 2);
    matrixToInvert(0,1) = ComplexNum(1, -1);
    matrixToInvert(1,1) = ComplexNum(-2, 2);
    matrixToInvert(2,1) = ComplexNum(2, -1);
    matrixToInvert(0,2) = ComplexNum(3, 0);
    matrixToInvert(1,2) = ComplexNum(5, 0);
    matrixToInvert(2,2) = ComplexNum(3, 1);

    Matrix inverse = inverseMatrix(&matrixToInvert);
    Matrix identity = matMul(&inverse, &matrixToInvert);
    std::cout << "inverse is\n" << inverse << std::endl;
    std::cout << "Result of matmul with matrix and its inverse is\n" << identity << std::endl;

    Matrix xData(6, 1);
    xData(0, 0) = 2;
    xData(1, 0) = 7;
    xData(2, 0) = 4;
    xData(3, 0) = 3;
    xData(4, 0) = 0;
    xData(5, 0) = -1;


    Matrix yData(6, 1);
    yData(0, 0) = 3;
    yData(1, 0) = 9;
    yData(2, 0) = 1;
    yData(3, 0) = 7;
    yData(4, 0) = -3;
    yData(5, 0) = -7;
    /*
    yData(0, 0) = 2;
    yData(1, 0) = 4;
    yData(2, 0) = 5;
    yData(3, 0) = 4;
    yData(4, 0) = 5;
    yData(5, 0) = 7;
    yData(6, 0) = 8;
    yData(7, 0) = 8;
    yData(8, 0) = 10;
    yData(9, 0) = 12;
     */


    LinearRegressor linreg1(&xData, &yData);
    std::cout << "Error is: " << linreg1.getError() << std::endl;
    Matrix testData(1, 1);
    testData(0, 0) = 825;
    ComplexNum prediction = linreg1.predict(&testData);
    std::cout << "Prediction is " << prediction << std::endl;
    Matrix regressionCoeff = linreg1.getRegressionCoefficients();
    std::cout << "Regression coefficients:\n" << regressionCoeff << std::endl;

    Matrix newMatrix = yData * ComplexNum(3, 0) ;
}
