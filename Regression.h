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

    // Setters


private:
    int numVars;
    double m_error;
    Matrix regressionCoefficients;
};


