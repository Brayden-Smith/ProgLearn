#include "Matrix.h"
#include <iostream>

Matrix::Matrix(int nr, int nc) : numRows(nr), numCols(nc), entryData(nr, std::vector<ComplexNum>(nc)) {
    if (nr <= 0 || nc <= 0) {
        throw std::invalid_argument("Matrix(int nr, int nc): Invalid size");
    }
}
Matrix::Matrix(const Matrix& matToCopy) : numRows(matToCopy.numRows), numCols(matToCopy.numCols), entryData(matToCopy.entryData) {}

Matrix::Matrix() : numRows(1), numCols(1) {}

// Arithmetic operators
Matrix Matrix::operator +(Matrix const& matrixToAdd) const {
    if (this->numRows != matrixToAdd.numRows || this->numCols != matrixToAdd.numCols) {
        throw std::invalid_argument("Matrix::operator +: Matrices not of same dimensions");
    }

    Matrix matrixToReturn(numRows, numCols);
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] + matrixToAdd.entryData[i][j];
        }
    }
    return matrixToReturn;
}

Matrix Matrix::operator*(ComplexNum const& scalar) {
    Matrix matrixToReturn(this->numRows, this->numCols);
    for (int i = 0; i < this->numRows; i++) {
        for (int j = 0; j < this->numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] * scalar;
        }
    }
    return matrixToReturn;
}

Matrix Matrix::operator*(double numToMul) {
    Matrix matrixToReturn(this->numRows, this->numRows);
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] * numToMul;
        }
    }
    return matrixToReturn;
}


Matrix Matrix::operator -(Matrix const& matrixToSub) const {
    if (this->numRows != matrixToSub.numRows || this->numCols != matrixToSub.numCols) {
        throw std::invalid_argument("Matrix::operator -: Matrices not of same dimensions");
    }

    Matrix matrixToReturn(matrixToSub.numRows, matrixToSub.numCols);
    for (int i = 0; i < this->numRows; i++) {
        for (int j = 0; j < this->numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] - matrixToSub.entryData[i][j];
        }
    }
    return matrixToReturn;
}

Matrix Matrix::operator -(ComplexNum const& numToSub) const {
    Matrix matrixToReturn(this->numRows, this->numCols);
    for (int i = 0; i < this->numRows; i++) {
        for (int j = 0; j < this->numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] - numToSub;
        }
    }
    return matrixToReturn;
}

// Equals
Matrix& Matrix::operator =(Matrix const& matrixToCopy) {
    /*
    if (this->numRows != matrixToCopy.numRows || this->numCols != matrixToCopy.numCols) {
        throw std::invalid_argument("Matrix::operator =: Matrices not of same dimensions");
    }
     */
    this->numRows = matrixToCopy.numRows;
    this->numCols = matrixToCopy.numCols;
    this->entryData = matrixToCopy.entryData;
    return *this;
}
bool Matrix::operator ==(const Matrix& otherMatrix) const {
    if (this->numRows != otherMatrix.numRows || this->numCols != otherMatrix.numCols) {
        return false;
    }

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            if (this->entryData[i][j] != otherMatrix.entryData[i][j]) {
                return false;
            }
        }
    }
    return true;
}

// Entry access operator
ComplexNum& Matrix::operator ()(int i, int j)  {
    if (this->numRows <= i || this->numCols <= j || i < 0 || j < 0) {
        throw std::invalid_argument("Matrix::operator (): Entry index out of bounds");
    }
    return this->entryData[i][j];
}

// Column accessor
Matrix Matrix::operator[](int n) const { // Honestly don't know which operator to use for this tried to make as intuitive as possible

    if (this->numCols <= n || n < 0) {
        throw std::invalid_argument("Matrix::operator[]: Index out of bounds");
    }

    Matrix result (this->numRows, 1);
    for(int i = 0; i < this->numRows; i++) {
        result(i, 0) = this->entryData[i][n];
    }
    return result;
}

//column mutator
void Matrix::columnAssign(int n, Matrix* M) {
    if(M->numCols > 1 || this->numRows < M->numRows) {
        throw std::invalid_argument("Matrix::columnAssign: Invalid column dimensions");
    }

    for(int i = 0; i < numRows; i++) {
        this->entryData[i][n] = (*M)(i,0);
    }
}

//Row accessor operator
Matrix Matrix::operator()(int n) const {
    if (this->numRows <= n || n < 0) {
        throw std::invalid_argument("Matrix::operator(): Index out of bounds");
    }

    Matrix result(1,this->numRows);
    for(int i = 0; i < this->numCols; i++) {
        result(0,i) = this->entryData[n][i];
    }
    return result;
}

//Row mutator
Matrix Matrix::rowAssign(int n, Matrix* M) {
    if(M->numRows > 1 || this->numCols < M->numCols) {
        throw std::invalid_argument("Matrix::rowAssign: Invalid dimensions ");
    }

    for(int i = 0; i < numCols; i++) {
        this->entryData[n][i] = (*M)(0,i);
    }
}

//Accessor methods
int Matrix::getNumRows() const{
    return this->numRows;
}
int Matrix::getNumCols() const {
    return this->numCols;
}


std::ostream &operator<<(std::ostream &outputStream, const Matrix& matrixToPrint) {

    for (int i = 0; i < matrixToPrint.numRows; i++) {
        outputStream << "| ";
        outputStream << matrixToPrint.entryData[i][0];
        for (int j  = 1; j < matrixToPrint.numCols - 1; j++) {
            outputStream << " " << matrixToPrint.entryData[i][j] << " ";
        }
        if (matrixToPrint.numCols > 1) {
            if (matrixToPrint.numCols == 2) {
                outputStream << ' ';
            }
            outputStream << matrixToPrint.entryData[i][matrixToPrint.numCols - 1];
        }
        outputStream << " |\n";
    }
    return outputStream;
}


Matrix operator*(ComplexNum const& scalar, Matrix const& matrix) {
    Matrix matrixToReturn(matrix.getNumRows(), matrix.getNumCols());
    for (int i = 0; i < matrix.getNumRows(); i++) {
        for (int j = 0; j < matrix.getNumCols(); j++) {
            matrixToReturn(i, j) = matrix.entryData[i][j] * scalar;
        }
    }
    return matrixToReturn;
}

Matrix operator*(double numToMul, Matrix const& matrix) {
    Matrix matrixToReturn(matrix.getNumRows(), matrix.getNumCols());
    for (int i = 0; i < matrix.getNumRows(); i++) {
        for (int j = 0; j < matrix.getNumCols(); j++) {
            matrixToReturn(i, j) = matrix.entryData[i][j] * numToMul;
        }
    }
    return matrixToReturn;
}

Matrix Matrix::conjugate() const {
    Matrix conjugate(this->numRows,this->numCols);
    for(int i = 0; i < this->numRows; i++) {

        for(int j = 0; j < this->numCols; j++) {
            conjugate(i,j) = this->entryData[i][j].getConjugate();
        }
    }
    return conjugate;
}


