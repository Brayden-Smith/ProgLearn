#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>



#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Waddress-of-temporary"








// BASIC MATRIX AND VECTOR FUNCTIONALITIES
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------




class ComplexNum {
public:
    // Constructor
    ComplexNum(double a = 0, double b = 0) : realPart(a), imagPart(b) {};

    // Equality operators
    ComplexNum& operator=(ComplexNum const& numToCopy) {
        if (this != &numToCopy) {
            this->realPart = numToCopy.realPart;
            this->imagPart = numToCopy.imagPart;
        }
        return *this;
    }
    ComplexNum& operator=(double numToCopy) {
        this->realPart = numToCopy;
        this->imagPart = 0;
        return *this;
    }
    ComplexNum& operator+=(ComplexNum const& numToAdd) {
        this->realPart += numToAdd.realPart;
        this->imagPart += numToAdd.imagPart;
        return *this;
    }

    bool operator==(const ComplexNum& otherNum) const {
        return (this->realPart == otherNum.realPart && this->imagPart == otherNum.imagPart);
    }
    bool operator==(double otherNum) const {
        return (this->realPart == otherNum && this->imagPart == 0);
    }
    bool operator!=(const ComplexNum& otherNum) const {
        return (!(this->realPart == otherNum.realPart && this->imagPart == otherNum.imagPart));
    }
    // Copy constructors

    ComplexNum(const ComplexNum& numToCopy) {
        this->realPart = numToCopy.realPart;
        this->imagPart = numToCopy.imagPart;
    }
    ComplexNum(double numToCopy) {
        this->realPart = numToCopy;
        this->imagPart = 0;
    }

    // Arithmetic operators

    ComplexNum operator+(ComplexNum const& numToAdd) const {
        return ComplexNum(this->realPart + numToAdd.realPart, this->imagPart + numToAdd.imagPart);
    }
    ComplexNum operator+(double numToAdd) const {
        return ComplexNum(this->realPart + numToAdd, this->imagPart);
    }

    ComplexNum operator*(ComplexNum const& numToMul) const {
        double a = this->realPart;
        double b = this->imagPart;
        double c = numToMul.realPart;
        double d = numToMul.imagPart;

        return ComplexNum((a * c) - (b * d), (a * d) + (b * c));
    }
    ComplexNum operator*(double numToMul) const {
        return ComplexNum (numToMul * this->realPart, numToMul * this->imagPart);
    }

    ComplexNum operator-(ComplexNum const& numToSub) const {
        return ComplexNum(this->realPart - numToSub.realPart, this->imagPart - numToSub.imagPart);
    }
    ComplexNum operator-(double numToSub) const {
        return ComplexNum(this->realPart - numToSub, this->imagPart);
    }

    ComplexNum operator/(double numToDiv) const {
        return ComplexNum(this->realPart / numToDiv, this->imagPart / numToDiv);
    }
    ComplexNum operator/(ComplexNum const& numToDiv) {
        double a2SquaredPlusb2Squared = (numToDiv.realPart * numToDiv.realPart) + (numToDiv.imagPart * numToDiv.imagPart);
        double a = (((this->realPart * numToDiv.realPart) + (this->imagPart * numToDiv.imagPart)) / a2SquaredPlusb2Squared);
        double b = ((this->imagPart * numToDiv.realPart) - (this->realPart * numToDiv.imagPart)) / a2SquaredPlusb2Squared;
        return ComplexNum(a,b);
    }


    friend ComplexNum operator+(double lhs, const ComplexNum& rhs);
    friend ComplexNum operator*(double lhs, const ComplexNum& rhs);

    // Method declarations
    double getMagnitude() const {
        return sqrt((realPart * realPart) + (imagPart * imagPart));
    };
    double getRealPart() const {
        return realPart;
    }
    double getImagPart() const {
        return imagPart;
    }
    ComplexNum getConjugate() const {
        ComplexNum numToReturn;
        numToReturn.realPart = this->realPart;
        numToReturn.imagPart = this->imagPart * -1;
        return numToReturn;
    }

    // Internal values
    double realPart;
    double imagPart;
};
// Friend functions for commutativity
ComplexNum operator+(double lhs, const ComplexNum& rhs) {
    return ComplexNum(lhs + rhs.realPart, rhs.imagPart);
}
ComplexNum operator*(double lhs, const ComplexNum& rhs) {
    return ComplexNum(lhs * rhs.realPart, lhs * rhs.imagPart);
}

// External complex conjugate methods added for easier use
ComplexNum complexConjugate(const ComplexNum& numToConjugate) {
    ComplexNum numToReturn;
    numToReturn.realPart = numToConjugate.realPart;
    numToReturn.imagPart = numToConjugate.imagPart * -1;
    return numToReturn;
}
ComplexNum complexConjugate(double numToConjugate) {
    ComplexNum numToReturn;
    numToReturn.realPart = numToConjugate;
    numToReturn.imagPart = 0;
    return numToReturn;
}
double magnitudeOfNumber(const ComplexNum& numToMag) {
    return sqrt((numToMag.realPart * numToMag.realPart) + (numToMag.imagPart * numToMag.imagPart));
}
double magnitudeOfNumber(double numToMag) {
    return abs(numToMag);
}
double complexNumToDouble(const ComplexNum& num) {
    if (num.imagPart != 0) {
        throw std::invalid_argument("Cannot convert complex number with non-zero imaginary part to double");
    }
    return num.realPart;
}



class Matrix {
public:
    // Constructors

    Matrix(int nr, int nc) : numRows(nr), numCols(nc), entryData(nr, std::vector<ComplexNum>(nc)) {
        if (nr <= 0 || nc <= 0) {
            throw std::invalid_argument("Invalid size");
        }
    }
    Matrix(const Matrix& matToCopy) : numRows(matToCopy.numRows), numCols(matToCopy.numCols), entryData(matToCopy.entryData) {}

    // Arithmetic operators
    Matrix operator +(Matrix const& matrixToAdd) const {
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
    Matrix operator *(ComplexNum scalar) {
        Matrix matrixToReturn(this->numRows, this->numCols);
        for (int i = 0; i < this->numRows; i++) {
            for (int j = 0; j < this->numCols; j++) {
                matrixToReturn.entryData[i][j] = this->entryData[i][j] * scalar;
            }
        }
        return matrixToReturn;
    }
    Matrix operator *(double numToMul) const {
        Matrix matrixToReturn(this->numRows, this->numRows);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                matrixToReturn.entryData[i][j] = this->entryData[i][j] * numToMul;
            }
        }
        return matrixToReturn;
    }
    friend Matrix operator*(double lhs, const Matrix& rhs);
    Matrix operator -(Matrix const& matrixToSub) const {
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
    Matrix& operator =(Matrix const& matrixToCopy) {
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
    bool operator ==(const Matrix& otherMatrix) const {
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
    ComplexNum& operator ()(int i, int j) { // Matrix access is one-indexed for more intuitive use. Feel free to change in a future commit if this proves annoying for QR decomposition
        if (this->numRows <= i || this->numCols <= j || i < 0 || j < 0) {
            throw std::invalid_argument("Entry index out of bounds");
        }
        return this->entryData[i][j];
    }

    // Column access operator
    Matrix operator[](int n) const { // Honestly don't know which operator to use for this tried to make as intuitive as possible

        if (this->numCols <= n || n < 0) {
            throw std::invalid_argument("Index out of bounds");
        }

        Matrix result (this->numRows, 1);
        for(int i = 0; i < this->numRows; i++) {
            result(i, 0) = this->entryData[i][n];
        }
        return result;
    }

    //column assign function
    void columnAssign(int n, std::vector<ComplexNum>* column) {
        if (this->numRows  != column->size()) {
            throw std::invalid_argument("Invalid column dimensions");
        }

        for(int i = 0; i < this->numRows; i++) {
            this->entryData[i][n] = (*column)[i];
        }
    }

    void columnAssign(int n, Matrix* M) {
        if(M->numCols > 1 || this->numRows < M->numRows) {
            throw std::invalid_argument("Invalid column dimensions");
        }

        for(int i = 0; i < numRows; i++) {
            this->entryData[i][n] = (*M)(i,0);
        }
    }


    // Internal data

    int numRows;
    int numCols;
    std::vector<std::vector<ComplexNum>> entryData; // Outer vector is row, inner vector is column. For example, entryData[1][2] accesses the second row, third column
};
Matrix operator*(double lhs, const Matrix &rhs) {
    Matrix matrixToReturn(rhs.numRows, rhs.numRows);
    for (int i = 0; i < rhs.numRows; i++) {
        for (int j = 0; i < rhs.numCols; j++) {
            matrixToReturn.entryData[i][j] = rhs.entryData[i][j] * lhs;
        }
    }
    return matrixToReturn;
}
Matrix matMul(Matrix const& lhs, Matrix const& rhs) {
    if (lhs.numCols != rhs.numRows) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }

    Matrix product(lhs.numRows, rhs.numCols);
    for (int i = 0; i < lhs.numRows; i++) {
        for (int j = 0; j < rhs.numCols; j++) {
            ComplexNum sum(0, 0);
            for (int k = 0; k < lhs.numCols; k++) {
                sum = sum + (lhs.entryData[i][k] * rhs.entryData[k][j]);
            }
            product.entryData[i][j] = sum;
        }
    }

    return product;
} // Slow, but functional

double realTraceOfMatrix(Matrix const& matrixToTrace) { // Only defined for real matrices at the moment
    if (matrixToTrace.numRows != matrixToTrace.numCols) {
        throw std::invalid_argument("Trace is only defined for square matrices");
    }

    double runningTrace = 0;
    for (int i = 0; i < matrixToTrace.numRows; i++) {
        double numberToAddToRunningTrace = complexNumToDouble(matrixToTrace.entryData[i][i]);
        runningTrace += numberToAddToRunningTrace;
    }

    return runningTrace;
}

Matrix conjTranspose(Matrix const& matrixToTranspose) {
    Matrix matrixToReturn(matrixToTranspose.numCols, matrixToTranspose.numRows);
    for (int i = 0; i < matrixToTranspose.numRows; i++) {
        for (int j = 0; j < matrixToTranspose.numCols; j++) {
            ComplexNum transposedEntry = complexConjugate(matrixToTranspose.entryData[i][j]);
            matrixToReturn.entryData[j][i] = transposedEntry;
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
    if (u.numCols != 1 || v.numCols != 1) {
        throw std::invalid_argument("Vectors must be nx1 matrices");
    }

    ComplexNum result(0, 0);
    for (int i = 0; i < u.numRows; i++) {
        ComplexNum vEntryConjugate = complexConjugate(v(i, 0));
        result += (u(i, 0) * vEntryConjugate);
    }
    return result;
}


// Method to normalize vectors in a matrix (alters matrix upon which it is called)
void normalizeVectorsInMatrix(Matrix* pointerToMatrix) {
    for (int i = 0; i < pointerToMatrix->numCols; i++) {

        double sumToNorm = 0;
        for (int j = 0; j < pointerToMatrix->numRows; j++) {
            double magnitude = magnitudeOfNumber(pointerToMatrix->entryData[j][i]);
            sumToNorm += (magnitude * magnitude);
        }

        double VecMagnitude = sqrt(sumToNorm);
        double oneOverMagnitude = 1/VecMagnitude;
        for (int j = 0; j < pointerToMatrix->numRows; j++) {
            pointerToMatrix->entryData[j][i] = pointerToMatrix->entryData[j][i] * oneOverMagnitude;
        }

    }
}

Matrix GramSchmidt(Matrix const& M) {
    Matrix result(M.numRows,M.numCols);

    for(int i = 0; i < M.numCols; i++) {
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







// MATRIX DECOMPOSITIONS AND EIGENVALUES
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------


// QR

// Eigenvalues

// SVD





int main() {
    // Complex number test asserts
    ComplexNum z(3, 2);
    assert(z.getConjugate().imagPart == -2 );
    assert ((2 + z).realPart == 5);
    assert((z * -7).imagPart == -14);
    assert((2 * z).realPart == 6);
    assert((z * 12).realPart == 36);
    assert((z * 0).realPart == 0);
    assert(z.getMagnitude() == sqrt(13));
    ComplexNum w(2, 3);
    assert((z*w).imagPart == 13);
    assert((z*w).realPart == 0);

    // Matrix test asserts
    Matrix A(2, 2);
    A(0,1) = ComplexNum(1, 0);
    assert(A(0, 1) == ComplexNum(1, 0));
    Matrix B(2, 2);
    B(0, 1) = ComplexNum(3, 4);
    Matrix C = A + B;
    assert(C(0, 1) == ComplexNum(4, 4));

    Matrix CTranspose = conjTranspose(C);
    assert(CTranspose(1, 0) == ComplexNum(4, -4));
    assert(CTranspose(0, 1) == ComplexNum(0, 0));

    Matrix matrixToNorm(3, 2);
    matrixToNorm(0, 0) = 3;
    matrixToNorm(1 ,0) = 4;
    matrixToNorm(2, 0) = 2;
    matrixToNorm(0, 1) = 1;
    matrixToNorm(1, 1) = ComplexNum(0, 2);
    matrixToNorm(2, 1) = 1;
    std::cout << frobeniusNorm(matrixToNorm) << std::endl;



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
    std::cout << "here";






};
#pragma clang diagnostic pop