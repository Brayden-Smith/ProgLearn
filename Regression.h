#include "LA.h"

Matrix ordinaryLeastSquaresCoefficients(Matrix* xData, Matrix* yData);


class LinearRegressor {
public:
    // Constructors
    LinearRegressor(Matrix* xData, Matrix* yData);

    // Operators
    LinearRegressor operator=(LinearRegressor const& rhs);

    // Methods
    ComplexNum predict(Matrix* xData);

    // Getters

    int getNumVars();
    double getError();
    Matrix getRegressionCoefficients();

    // Setters


private:
    int numVars;
    double m_error;
    Matrix regressionCoefficients;
};

ComplexNum sigmoid (ComplexNum numberToSigmoid);

void yDataProbTransform (Matrix* yDataToTransform);

ComplexNum logit (ComplexNum prob);

class OlsLogisticRegressor {
public:

    OlsLogisticRegressor();
    void fit(Matrix* xData, Matrix* yData);
    Matrix predict(Matrix* xDataToPred);
    Matrix getRegressionCoefficients();


private:
    Matrix regressionCoefficients;
    int numVars;
    bool isFit;
};

class RidgeRegressor {
public:


private:

};
