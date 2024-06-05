#include "Matrix.h"

Matrix::Matrix(int nr, int nc) : numRows(nr), numCols(nc), entryData(nr, std::vector<ComplexNum>(nc)) {
    if (nr <= 0 || nc <= 0) {
        throw std::invalid_argument("Invalid size");
    }
}
Matrix::Matrix(const Matrix& matToCopy) : numRows(matToCopy.numRows), numCols(matToCopy.numCols), entryData(matToCopy.entryData) {}

// Arithmetic operators
Matrix Matrix::operator +(Matrix const& matrixToAdd) const {
    if (this->numRows != matrixToAdd.numRows || this->numCols != matrixToAdd.numCols) {
        throw std::invalid_argument("Matrices not of same dimensions");
    }

    Matrix matrixToReturn(numRows, numCols);
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] + matrixToAdd.entryData[i][j];
        }
    }
    return matrixToReturn;
}
Matrix Matrix::operator *(ComplexNum const& scalar) {
    Matrix matrixToReturn(this->numRows, this->numCols);
    for (int i = 0; i < this->numRows; i++) {
        for (int j = 0; j < this->numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] * scalar;
        }
    }
    return matrixToReturn;
}
Matrix Matrix::operator *(double numToMul) const {
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
        throw std::invalid_argument("Matrices not of same dimensions");
    }

    Matrix matrixToReturn(matrixToSub.numRows, matrixToSub.numCols);
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            matrixToReturn.entryData[i][j] = this->entryData[i][j] - matrixToSub.entryData[i][j];
        }
    }
    return matrixToReturn;
}

// Equals
Matrix& Matrix::operator =(Matrix const& matrixToCopy) {
    if (this->numRows != matrixToCopy.numRows || this->numCols != matrixToCopy.numCols) {
        throw std::invalid_argument("Matrices not of same dimensions");
    }

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            this->entryData[i][j] = matrixToCopy.entryData[i][j];
        }
    }
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
        throw std::invalid_argument("Entry index out of bounds");
    }
    return this->entryData[i][j];
}

// Column accessor
Matrix Matrix::operator[](int n) const { // Honestly don't know which operator to use for this tried to make as intuitive as possible

    if (this->numCols <= n || n < 0) {
        throw std::invalid_argument("Index out of bounds");
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
        throw std::invalid_argument("Invalid column dimensions");
    }

    for(int i = 0; i < numRows; i++) {
        this->entryData[i][n] = (*M)(i,0);
    }
}

//Row accessor operator
Matrix Matrix::operator()(int n) const {
    if (this->numRows <= n || n < 0) {
        throw std::invalid_argument("Index out of bounds");
    }

    Matrix result(1,this->numRows);
    for(int i = 0; i < this->numCols; i++) {
        result(0,i) = this->entryData[n][i];
    }
}

//Row mutator
Matrix Matrix::rowAssign(int n, Matrix* M) {
    if(M->numRows > 1 || this->numCols < M->numCols) {
        throw std::invalid_argument("Invalid dimensions");
    }

    for(int i = 0; i < numCols; i++) {
        this->entryData[n][i] = (*M)(0,i);
    }
}

//Accessor methods
double Matrix::getNumRows() const{
    return this->numRows;
}
double Matrix::getNumCols() const {
    return this->numCols;
}

Matrix operator*(double lhs, const Matrix &rhs) {
    Matrix matrixToReturn(rhs.numRows, rhs.numRows);
    for (int i = 0; i < rhs.numRows; i++) {
        for (int j = 0; i < rhs.numCols; j++) {
            matrixToReturn.entryData[i][j] = rhs.entryData[i][j] * lhs;
        }
    }
    return matrixToReturn;
}
Matrix matMul(Matrix* lhs, Matrix* rhs) {
    if (lhs->getNumCols() != rhs->getNumRows()) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }

    Matrix product(lhs->getNumRows(), rhs->getNumCols());
    for (int i = 0; i < lhs->getNumRows(); i++) {
        for (int j = 0; j < rhs->getNumCols(); j++) {
            ComplexNum sum(0, 0);
            for (int k = 0; k < lhs->getNumCols(); k++) {
                sum = sum + ((*lhs)(i,k) * (*rhs)(k,j));
            }
            product(i,j) = sum;
        }
    }

    return product;
}

double realTraceOfMatrix(Matrix* matrixToTrace) { // Only defined for real matrices at the moment
    if (matrixToTrace->getNumRows() != matrixToTrace->getNumCols()) {
        throw std::invalid_argument("Trace is only defined for square matrices");
    }

    double runningTrace = 0;
    for (int i = 0; i < matrixToTrace->getNumRows(); i++) {
        double numberToAddToRunningTrace = complexNumToDouble((*matrixToTrace)(i,i));
        runningTrace += numberToAddToRunningTrace;
    }

    return runningTrace;
}

Matrix conjTranspose(Matrix* matrixToTranspose) {
    Matrix matrixToReturn(matrixToTranspose->getNumCols(), matrixToTranspose->getNumRows());
    for (int i = 0; i < matrixToTranspose->getNumRows(); i++) {
        for (int j = 0; j < matrixToTranspose->getNumCols(); j++) {
            ComplexNum transposedEntry = complexConjugate((*matrixToTranspose)(i,j));
            matrixToReturn(j,i) = transposedEntry;
        }
    }
    return matrixToReturn;
}

double frobeniusNorm(Matrix const& matrixToNorm) {
    Matrix conjugateTransposeMatrix = conjTranspose(matrixToNorm);
    Matrix matrixStarMatrix = matMul(conjugateTransposeMatrix, matrixToNorm);
    //std::cout << "Matrix mul has dimensions " << matrixStarMatrix.numRows << "x" << matrixStarMatrix.numCols << std::endl;
    double trace = realTraceOfMatrix(matrixStarMatrix);
    return sqrt(trace);
}

ComplexNum InnerProduct(Matrix& u, Matrix& v) {
    if (u.getNumCols() != 1 || v.getNumCols() != 1) {
        throw std::invalid_argument("Vectors must be nx1 matrices");
    }

    ComplexNum result(0, 0);
    for (int i = 0; i < u.getNumRows(); i++) {
        ComplexNum vEntryConjugate = complexConjugate(v(i, 0));
        result += (u(i, 0) * vEntryConjugate);
    }
    return result;
}


// Method to normalize vectors in a matrix (alters matrix upon which it is called)
void normalizeVectorsInMatrix(Matrix* pointerToMatrix) {
    for (int i = 0; i < pointerToMatrix->getNumCols(); i++) {

        double Norm = 0;
        for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
            double magnitude = magnitudeOfNumber((*pointerToMatrix)(j,i));
            Norm += (magnitude * magnitude);
        }

        Norm = sqrt(Norm);
        double oneOverMagnitude = 1/Norm;
        for (int j = 0; j < pointerToMatrix->getNumRows(); j++) {
            (*pointerToMatrix)(j,i) = (*pointerToMatrix)(j,i) * oneOverMagnitude;
        }

    }
}

Matrix GramSchmidt(Matrix const& M) {
    Matrix result(M.getNumRows(),M.getNumCols());

    for(int i = 0; i < M.getNumCols(); i++) {
        Matrix Vk = M[i];
        Matrix Uk = Vk;
        for(int j = 0; j < i; j++) {

            Matrix Un = result[j];
            ComplexNum VdotU = InnerProduct(Vk,Un);
            ComplexNum UdotU = InnerProduct(Un,Un);

            //compute the projection
            ComplexNum quotient = (VdotU / UdotU);
            Un = Un * quotient;
            Uk = Uk - Un;
        }
        result.columnAssign(i,&Uk);
    }

    normalizeVectorsInMatrix(&result);
    return result;
}