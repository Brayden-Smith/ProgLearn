#include "Regression.h"

Matrix ordinaryLeastSquaresCoefficients(Matrix* xData, Matrix* yData) { // XData must be invertible. Entries in original XData are horizontal (one row is one observation)
// Internally, one column is one observation, but not until after the transpose
    if (xData->getNumRows() <= xData->getNumCols()) {
        throw std::invalid_argument("ordinaryLeastSquaresCoefficients: Must have more observations than predictors"); // To decrease the chance of XX^T being singular
    }

    if (yData->getNumCols() > 1) {
        throw std::invalid_argument("ordinaryLeastSquaresCoefficients: Invalid y vector!");
    }

    Matrix columnOfOnes(xData->getNumRows(), 1);
    for (int i = 0; i < xData->getNumRows(); i++) {
        columnOfOnes(i, 0) = 1;
    }
    Matrix xDataA = createAugmentedMatrix(&columnOfOnes, xData);
    Matrix xDataACorrected = transpose(&xDataA);
    Matrix* xDataAPointer = &xDataACorrected;
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
    if (MSE < 1e-9) {
        m_error = 0;
    }
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

Matrix LinearRegressor::getRegressionCoefficients() {
    return this->regressionCoefficients;
}


// Logistic regression

ComplexNum sigmoid (ComplexNum numberToSigmoid) {
    if (numberToSigmoid.getImagPart() != 0) {
        throw std::invalid_argument("sigmoid: Sigmoids are only defined for real numbers!"); // This can change later, but at the moment logistic regression will only be defined for real numbers
    }

    double sigmoid = 1/(1 + std::exp(numberToSigmoid.getRealPart() * -1));

    return {sigmoid, 0};
}

void yDataProbTransform(Matrix* yDataToTransform) {
    for (int i = 0; i < yDataToTransform->getNumRows(); i++) {
        if ((*yDataToTransform)(i, 0) == 0) {
            (*yDataToTransform)(i, 0) = ComplexNum(0.0001, 0);
        } else if ((*yDataToTransform)(i, 0) == 1) {
            (*yDataToTransform)(i, 0) = ComplexNum(0.9999, 0);
        }
    }
}

ComplexNum logit(ComplexNum prob) {
    if (prob.getImagPart() != 0) {
        throw std::invalid_argument("logit: Logits are only defined for real numbers!"); // This can change later, but at the moment logistic regression will only be defined for real numbers
    }

    ComplexNum odds = prob/((-1 * prob) + 1);
    double logOdds = std::log(odds.getRealPart());
    ComplexNum cnumToRet(logOdds, 0);
    return cnumToRet;
}

OlsLogisticRegressor::OlsLogisticRegressor() : regressionCoefficients(), numVars(0), isFit(false) {}

void OlsLogisticRegressor::fit(Matrix* xData, Matrix* yData) {
    numVars = xData->getNumCols();
    Matrix copyOfYData = (*yData);
    yDataProbTransform(&copyOfYData);

    Matrix yLogits(yData->getNumRows(), 1);

    for (int i = 0; i < yData->getNumRows(); i++) {
        ComplexNum zToTransform = copyOfYData(i, 0);
        yLogits(i, 0) = logit(zToTransform);
    }

    Matrix coefficients = ordinaryLeastSquaresCoefficients(xData, &yLogits);

    isFit = true;
    regressionCoefficients = coefficients;

}

Matrix OlsLogisticRegressor::predict(Matrix *xDataToPred) {
    Matrix yHat(xDataToPred->getNumRows(), 1);

    for (int i = 0; i < xDataToPred->getNumRows(); i++) {
        Matrix xRow(1, xDataToPred->getNumCols());
        for (int j = 0; j < xDataToPred->getNumCols(); j++) {
            xRow(0, j) = (*xDataToPred)(i, j);
        }

        Matrix columnOfOnes(1, 1);
        columnOfOnes(0, 0) = 1;

        Matrix xDataA = createAugmentedMatrix(&columnOfOnes, &xRow);
        //std::cout << "Predict augmented matrix is\n" << xDataA << std::endl;
        Matrix betaMinizerTimesData = matMul(&xDataA, &regressionCoefficients);
        ComplexNum numToSigmoid = betaMinizerTimesData(0, 0);

        ComplexNum sigmoidRes = sigmoid(numToSigmoid);

        yHat(i, 0) = sigmoidRes;
    }

    return yHat;
}

Matrix OlsLogisticRegressor::getRegressionCoefficients() {
    return regressionCoefficients;
}



