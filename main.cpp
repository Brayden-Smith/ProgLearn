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

    Matrix Fill(2,1);
    Fill(0,0) = 2;
    Fill(1,0) = 3;
    //auto t = frontFillVec(Fill,5,ComplexNum(5,5));
    //std::cout << t;

    Matrix matrixToDecomp(4, 3);
    matrixToDecomp(0,0) = ComplexNum(1, 0);
    matrixToDecomp(1,0) = ComplexNum(1, 0);
    matrixToDecomp(2,0) = ComplexNum(1, 0);
    matrixToDecomp(3,0) = ComplexNum(1, 0);
    matrixToDecomp(0,1) = ComplexNum(-1, 0);
    matrixToDecomp(1,1) = ComplexNum(4, 0);
    matrixToDecomp(2,1) = ComplexNum(4, 0);
    matrixToDecomp(3,1) = ComplexNum(-1, 0);
    matrixToDecomp(0,2) = ComplexNum(4, 0);
    matrixToDecomp(1,2) = ComplexNum(-2, 0);
    matrixToDecomp(2,2) = ComplexNum(2, 0);
    matrixToDecomp(3,2) = ComplexNum(0, 0);
    //std::cout << "the Q of QR Decomposition\n" << QRDecomp(matrixToDecomp)[0];
    //std::cout << "Matrix to decomp\n" << matrixToDecomp;
    //std::cout << "The R of QR decomp\n" << QRDecomp(matrixToDecomp)[1];


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


    // OLS regression test
    LinearRegressor linreg1(&xData, &yData);
    //std::cout << "Error is: " << linreg1.getError() << std::endl;
    Matrix testData(1, 1);
    testData(0, 0) = 825;
    ComplexNum prediction = linreg1.predict(&testData);
    //std::cout << "Prediction is " << prediction << std::endl;
    Matrix regressionCoeff = linreg1.getRegressionCoefficients();
    //std::cout << "Regression coefficients:\n" << regressionCoeff << std::endl;

    // Eigenvalues test

    Matrix wantEigenvalues(3, 3);

    /*
    wantEigenvalues(0, 0) = 6;
    wantEigenvalues(1, 0) = 2;
    wantEigenvalues(2, 0) = 1;
    wantEigenvalues(0, 1) = 2;
    wantEigenvalues(1, 1) = 3;
    wantEigenvalues(2, 1) = 1;
    wantEigenvalues(0, 2) = 1;
    wantEigenvalues(1, 2) = 1;
    wantEigenvalues(2, 2) = 1;
     */

    wantEigenvalues(0, 0) = 1;
    wantEigenvalues(1, 0) = 3;
    wantEigenvalues(2, 0) = 2;
    wantEigenvalues(0, 1) = 2;
    wantEigenvalues(1, 1) = 2;
    wantEigenvalues(2, 1) = 1;
    wantEigenvalues(0, 2) = 3;
    wantEigenvalues(1, 2) = 1;
    wantEigenvalues(2, 2) = 3;

    //std::cout << "The R of QR decomp\n" << QRDecomp(wantEigenvalues)[1];


    std::vector<ComplexNum> eigens = eigenvalues(&wantEigenvalues);
    std::cout << "Eigenvalues" << std::endl;
    std::cout << "eigens length is " << eigens.size() << std::endl;
    for (int i = 0; i < eigens.size(); i++) {
        std::cout << eigens[i] << std::endl;
    }

    //expected value test
    Matrix E(5,1);
    E(0,0) = ComplexNum(5,0);
    E(1,0) = ComplexNum(2,0);
    E(2,0) = ComplexNum(3,0);
    E(3,0) = ComplexNum(2,0);
    E(4,0) = ComplexNum(1,0);
    std::cout<< "\nExpected value: " << expectedValue(&E) << std::endl;

    //covariance matrix test
    /*
    Matrix covariance(4,3);
    covariance(0,0) = ComplexNum(1, 0);
    covariance(1,0) = ComplexNum(1, 0);
    covariance(2,0) = ComplexNum(1, 0);
    covariance(3,0) = ComplexNum(1, 0);
    covariance(0,1) = ComplexNum(-1, 0);
    covariance(1,1) = ComplexNum(4, 0);
    covariance(2,1) = ComplexNum(4, 0);
    covariance(3,1) = ComplexNum(-1, 0);
    covariance(0,2) = ComplexNum(4, 0);
    covariance(1,2) = ComplexNum(-2, 0);
    covariance(2,2) = ComplexNum(2, 0);
    covariance(3,2) = ComplexNum(0, 0);
    covariance = covarianceMatrix(&covariance);
    std::cout << covariance;
     */

}
