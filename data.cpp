#include "data.h"
#include <random>

Matrix shuffleRowOrder(Matrix& matrix) {

}

std::vector<Matrix> trainTestSplit(Matrix& data, double percent) {
    if (percent > 1 || percent < 0) {
        throw std::invalid_argument("trainTestSplit: percent must be between 0 and 1!");
    }

    int numRows = data.getNumRows();
    int numTrainRows = std::floor(percent * numRows);
    int numTestRows = numRows - numTrainRows;

    Matrix trainMatrix(numTrainRows, data.getNumCols());
    Matrix testMatrix(numTestRows, data.getNumCols());


    for (int i = 0; i < numTrainRows; i++) {
        for (int j = 0; j < data.getNumCols(); j++) {
            trainMatrix(i, j) = data(i, j);
        }
    }

    for (int i = numTrainRows; i < numRows; i++) {
        for (int j = 0; j < data.getNumCols(); j++) {
            testMatrix(i - numTrainRows, j) = data(i, j);
        }
    }

    return {trainMatrix, testMatrix};
}
