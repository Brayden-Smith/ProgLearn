#include "Regression.h"




Matrix ordinaryLeastSquaresCoefficients(Matrix* xData, Matrix* yData) { // XData must be invertible. Entries in original XData are horizontal (one row is one observation)
// Internally, one column is one observation, but not until after the transpose
    if (xData->getNumRows() <= xData->getNumCols()) {
        throw std::invalid_argument("Must have more observations than predictors"); // To decrease the chance of XX^T being singular
    }

    if (yData->getNumCols() > 1) {
        throw std::invalid_argument("Invalid y vector!");
    }

    Matrix columnOfOnes(xData->getNumRows(), 1);
    for (int i = 0; i < xData->getNumRows(); i++) {
        columnOfOnes(i, 0) = 1;
    }
    Matrix xDataA = createAugmentedMatrix(&columnOfOnes, xData);
    Matrix xDataAtrans = transpose(&xDataA);
    Matrix* xDataAPointer = &xDataAtrans;

    Matrix xDataTranspose = conjTranspose(xDataAPointer);
    Matrix xDataXDataTranspose = matMul(xDataAPointer, &xDataTranspose);
    Matrix xDataXDataTransposeInverse = inverseMatrix(&xDataXDataTranspose);
    Matrix xDataY = matMul(xDataAPointer, yData);
    Matrix betaMinimizer = matMul(&xDataXDataTransposeInverse, &xDataY);
    return betaMinimizer;
}

LinearRegressor::LinearRegressor(Matrix* xData, Matrix* yData) : numVars(xData->getNumRows()), regressionCoefficients(ordinaryLeastSquaresCoefficients(xData, yData)) {
    Matrix columnOfOnes(xData->getNumRows(), 1);
    for (int i = 0; i < xData->getNumRows(); i++) {
        columnOfOnes(i, 0) = 1;
    }
    Matrix xDataA = createAugmentedMatrix(&columnOfOnes, xData);

    //std::cout << "yData:\n" << *yData << std::endl;
    //std::cout << "xDataA\n" << xDataA << std::endl;
    //std::cout << "Regression coefficients:\n" << regressionCoefficients << std::endl;
    Matrix yDataHat = matMul(&xDataA, &regressionCoefficients);

    double MSE = 0;

    //std::cout << "Y hat is:\n" << yDataHat << std::endl;
    for (int i = 0; i < yData->getNumRows(); i++) {
        MSE += (magnitudeOfNumber((*yData)(i, 0) - yDataHat(i, 0))) * (magnitudeOfNumber((*yData)(i, 0) - yDataHat(i, 0)));
    }

    MSE = MSE/yData->getNumRows();
    m_error = MSE;
}

ComplexNum LinearRegressor::predict(Matrix* xData) {
    Matrix columnOfOnes(xData->getNumRows(), 1);
    for (int i = 0; i < xData->getNumRows(); i++) {
        columnOfOnes(i, 0) = 1;
    }

    Matrix xDataA = createAugmentedMatrix(&columnOfOnes, xData);
    //std::cout << "Predict augmented matrix is\n" << xDataA << std::endl;
    Matrix betaMinizerTimesData = matMul(&xDataA, &regressionCoefficients);
    return betaMinizerTimesData(0, 0);
}

double LinearRegressor::getError() {
    return this->m_error;
}