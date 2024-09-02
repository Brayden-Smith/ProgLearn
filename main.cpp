#include <iostream>
#include "Matrix.h"
#include <assert.h>
#include "LA.h"
#include "ComplexNum.h"
#include "Regression.h"
#include <opencv2/opencv.hpp>
#include <windows.h>
#include "data.h"
#include "eigenfaces.h"
#include <filesystem>

int main() {





    // Complex number test
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

    //std::cout << ComplexNum(-4, -4) << std::endl;

    Matrix Fill(2,1);
    Fill(0,0) = 2;
    Fill(1,0) = 3;
    //auto t = frontFillVec(Fill,5,ComplexNum(5,5));
    //std::cout << t;

    /*
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
    */


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



    // OLS regression test
    LinearRegressor linreg1(&xData, &yData);
    //std::cout << "Error is: " << linreg1.getError() << std::endl;
    Matrix testData(1, 1);
    testData(0, 0) = 825;
    ComplexNum prediction = linreg1.predict(&testData);
    //@std::cout << "Prediction is " << prediction << std::endl;
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
    /*
    wantEigenvalues(0, 0) = 1;
    wantEigenvalues(1, 0) = 3;
    wantEigenvalues(2, 0) = 2;
    wantEigenvalues(0, 1) = 2;
    wantEigenvalues(1, 1) = 2;
    wantEigenvalues(2, 1) = 1;
    wantEigenvalues(0, 2) = 3;
    wantEigenvalues(1, 2) = 1;
    wantEigenvalues(2, 2) = 3;
    */
    wantEigenvalues(0, 0) = 1;
    wantEigenvalues(1, 0) = 2;
    wantEigenvalues(2, 0) = 3;
    wantEigenvalues(0, 1) = 2;
    wantEigenvalues(1, 1) = 4;
    wantEigenvalues(2, 1) = 5;
    wantEigenvalues(0, 2) = 3;
    wantEigenvalues(1, 2) = 5;
    wantEigenvalues(2, 2) = 6;

    Matrix matToSVD(2, 3);
    matToSVD(0, 0) = 1;
    matToSVD(1, 0) = 4;
    matToSVD(0, 1) = 2;
    matToSVD(1, 1) = 5;
    matToSVD(0, 2) = 3;
    matToSVD(1, 2) = 6;

    matToSVD = transpose(&matToSVD);

    //std::cout << "The R of QR decomp\n" << QRDecomp(wantEigenvalues)[1];


    //std::vector<ComplexNum> eigens = eigenvalues(&wantEigenvalues);
    //std::cout << "Eigenvalues" << std::endl;
    //std::cout << "eigens length is " << eigens.size() << std::endl;


    //expected value test
    Matrix E(1,5);
    E(0,0) = ComplexNum(5,0);
    E(0,1) = ComplexNum(2,0);
    E(0,2) = ComplexNum(3,0);
    E(0,3) = ComplexNum(2,0);
    E(0,4) = ComplexNum(1,0);
    assert(expectedValue(&E) == ComplexNum(2.6,0));


    //covariance test
    Matrix x2Data(4,1);
    x2Data(0,0) = ComplexNum(1,0);
    x2Data(1,0) = ComplexNum(-1,0);
    x2Data(2,0) = ComplexNum(4,0);
    x2Data(3,0) = ComplexNum(0,0);
    Matrix y2Data(4,1);
    y2Data(0,0) = ComplexNum(1,0);
    y2Data(1,0) = ComplexNum(-1,0);
    y2Data(2,0) = ComplexNum(4,0);
    y2Data(3,0) = ComplexNum(0,0);
    //std::cout<<covariance(&x2Data,&y2Data);



    //covariance matrix test

    Matrix covariance(4,3);
    covariance(0,0) = ComplexNum(-3, 0);
    covariance(1,0) = ComplexNum(1, 0);
    covariance(2,0) = ComplexNum(1, 0);
    covariance(3,0) = ComplexNum(-4, 0);
    covariance(0,1) = ComplexNum(-3, 0);
    covariance(1,1) = ComplexNum(-1, 0);
    covariance(2,1) = ComplexNum(8, 0);
    covariance(3,1) = ComplexNum(8, 0);
    covariance(0,2) = ComplexNum(-2, 0);
    covariance(1,2) = ComplexNum(1, 0);
    covariance(2,2) = ComplexNum(3, 0);
    covariance(3,2) = ComplexNum(-3, 0);
    covariance = covarianceMatrix(&covariance);
    std::cout << covariance;




    // Logistic regression tests
    //std::cout << "LOG REG\n" << std::endl;
    Matrix logregx(7, 2);
    logregx(0, 0) = 22;
    logregx(1, 0) = 25;
    logregx(2, 0) = 47;
    logregx(3, 0) = 52;
    logregx(4, 0) = 46;
    logregx(5, 0) = 56;
    logregx(6, 0) = 48;

    logregx(0, 1) = 20000;
    logregx(1, 1) = 35000;
    logregx(2, 1) = 50000;
    logregx(3, 1) = 45000;
    logregx(4, 1) = 30000;
    logregx(5, 1) = 60000;
    logregx(6, 1) = 70000;

    Matrix logregy(7, 1);
    logregy(0, 0) = 0;
    logregy(1, 0) = 0;
    logregy(2, 0) = 1;
    logregy(3, 0) = 1;
    logregy(4, 0 ) = 0;
    logregy(5, 0) = 1;
    logregy(6, 0) = 1;

    LinearRegressor linreg2(&logregx, &logregy);
    //std::cout << "Ling reg coefficients for this dataset are:\n" << linreg2.getRegressionCoefficients() << std::endl;

    OlsLogisticRegressor logreg;
    logreg.fit(&logregx, &logregy);

    Matrix xpred(2, 2);
    xpred(0, 0) = 45;
    xpred(0, 1) = 65000;
    xpred(1, 0) = 23;
    xpred(1, 1) = 60000;

    Matrix ypred = logreg.predict(&xpred);
    //std::cout << "ypred is:\n" << ypred << std::endl;
    //std::cout << "Regression coefficients are:\n" << logreg.getRegressionCoefficients() << std::endl;

    //std::cout << "Wanteigenvalues:\n" << wantEigenvalues << std::endl;

    std::vector<ComplexNum> core;
    std::vector<Matrix> decomp = singularValueDecomp(&matToSVD);

    std::cout << "Matrix to SVD:\n" << matToSVD << std::endl;

    std::cout << "U:\n" << decomp[0] << std::endl;
    std::cout << "Sigma:\n" << decomp[1] << std::endl;
    std::cout << "V T:\n" << decomp[2] << std::endl;

    Matrix intermed = matMul(&decomp[1], &decomp[2]);
    Matrix originalMatrix = matMul(&decomp[0], &intermed);
    std::cout << "Orig:\n" << originalMatrix << std::endl;

    cv::Mat testMat = cv::Mat::zeros(3, 3, CV_8U);
    std::cout << "Created OpenCV Mat: " << testMat << std::endl;


    std::string path = "C:\\Users\\chris\\Desktop\\Faces\\one.pgm";
    Matrix image = importGrayscaleImage(path);

    std::string hi = "Hello";
    displayImage(hi, image);


}
