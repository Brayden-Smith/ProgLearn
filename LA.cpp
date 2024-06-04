#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>


// Classes

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Waddress-of-temporary"
class ComplexNum {
public:
    // Constructor
    ComplexNum(double a = 0, double b = 0) : realPart(a), complexPart(b) {};

    // Equality operators
    ComplexNum& operator=(ComplexNum const& numToCopy) {
        if (this != &numToCopy) {
            this->realPart = numToCopy.realPart;
            this->complexPart = numToCopy.complexPart;
        }
        return *this;
    }
    ComplexNum& operator=(double numToCopy) {
        this->realPart = numToCopy;
        this->complexPart = 0;
        return *this;
    }
    ComplexNum& operator+=(ComplexNum const& numToAdd) {
        this->realPart += numToAdd.realPart;
        this->complexPart += numToAdd.complexPart;
        return *this;
    }

    bool operator==(const ComplexNum& otherNum) const {
        return (this->realPart == otherNum.realPart && this->complexPart == otherNum.complexPart);
    }
    bool operator!=(const ComplexNum& otherNum) const {
        return (!(this->realPart == otherNum.realPart && this->complexPart == otherNum.complexPart));
    }
    // Copy constructors

    ComplexNum(const ComplexNum& numToCopy) {
        this->realPart = numToCopy.realPart;
        this->complexPart = numToCopy.complexPart;
    }
    ComplexNum(double numToCopy) {
        this->realPart = numToCopy;
        this->complexPart = 0;
    }

    // Arithmetic operators

    ComplexNum operator+(ComplexNum const& numToAdd) const {
        return ComplexNum(this->realPart + numToAdd.realPart, this->complexPart + numToAdd.complexPart);
    }
    ComplexNum operator+(double numToAdd) const {
        return ComplexNum(this->realPart + numToAdd, this->complexPart);
    }

    ComplexNum operator*(ComplexNum const& numToMul) const {
        double a = this->realPart;
        double b = this->complexPart;
        double c = numToMul.realPart;
        double d = numToMul.complexPart;

        return ComplexNum((a * c) - (b * d), (a * d) + (b * c));
    }
    ComplexNum operator*(double numToMul) const {
        return ComplexNum (numToMul * this->realPart, numToMul * this->complexPart);
    }

    ComplexNum operator-(ComplexNum const& numToSub) const {
        return ComplexNum(this->realPart - numToSub.realPart, this->complexPart - numToSub.complexPart);
    }
    ComplexNum operator-(double numToSub) const {
        return ComplexNum(this->realPart - numToSub, this->complexPart);
    }

    ComplexNum operator/(double numToDiv) const {
        return ComplexNum(this->realPart / numToDiv, this->complexPart / numToDiv);
    }
    ComplexNum operator/(ComplexNum const& numToDiv) {
        double a2SquaredPlusb2Squared = (numToDiv.realPart * numToDiv.realPart) + (numToDiv.complexPart * numToDiv.complexPart);
        double a = (((this->realPart * numToDiv.realPart) + (this->complexPart * numToDiv.complexPart)) / a2SquaredPlusb2Squared);
        double b = ((this->complexPart * numToDiv.realPart) - (this->realPart * numToDiv.complexPart)) / a2SquaredPlusb2Squared;
        return ComplexNum(a,b);
    }
    // todo(?) division of two complex numbers


    friend ComplexNum operator+(double lhs, const ComplexNum& rhs);
    friend ComplexNum operator*(double lhs, const ComplexNum& rhs);

    // Method declarations
    double getMagnitude() const {
        return sqrt((realPart * realPart) + (complexPart * complexPart));
    };
    double getRealPart() const {
        return realPart;
    }
    double getComplexPart() const {
        return complexPart;
    }
    ComplexNum getConjugate() const {
        ComplexNum numToReturn;
        numToReturn.realPart = this->realPart;
        numToReturn.complexPart = this->complexPart * -1;
        return numToReturn;
    }

    // Internal values
    double realPart;
    double complexPart;
};
// Friend functions for commutativity
ComplexNum operator+(double lhs, const ComplexNum& rhs) {
    return ComplexNum(lhs + rhs.realPart, rhs.complexPart);
}
ComplexNum operator*(double lhs, const ComplexNum& rhs) {
    return ComplexNum(lhs * rhs.realPart, lhs * rhs.complexPart);
}



class Matrix {
public:
    // Constructors

    Matrix(int nr, int nc) : numRows(nr), numCols(nc), entryData(nr, std::vector<ComplexNum>(nc)) {
        if (nr <= 0 || nc <= 0) {
            throw std::invalid_argument("Invalid size");
        }
    }
    Matrix(const Matrix& matToCopy) : numRows(matToCopy.numRows), numCols(matToCopy.numCols), entryData(matToCopy.entryData) {

    }

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
    Matrix operator *(double numToMul) const {
        Matrix matrixToReturn(this->numRows, this->numRows);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; i < numCols; j++) {
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
            for (int j = 0; i < numCols; j++) {
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
            for (int j = 0; i < numCols; j++) {
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
        if (this->numRows < i - 1 || this->numCols < j - 1 || i < 1 || j < 1) {
            throw std::invalid_argument("Entry index out of bounds");
        }
        return this->entryData[i - 1][j - 1];
    }

    //column access operator
    std::vector<ComplexNum> operator[](int n) { //Honestly don't know which operator to use for this tried to make as intuitive as possible

        if (this->numCols < n || n < 1) {
            throw std::invalid_argument("Index out of bounds");
        }

        std::vector<ComplexNum> result;

        for(int i = 1; i <= this->numRows; i++) {
            result.push_back(this->operator()(i,n));
        }
        return result;
    }

    //column assign function
    void columnAssign(int n, std::vector<ComplexNum>* column) {
        if (this->numRows  != column->size()) {
            throw std::invalid_argument("Invalid dimensions");
        }


        for(int i = 0; i < numRows; i++) {
            this->entryData[i][n-1] = (*column)[i];
        }
    }


    // Internal data

    int numRows;
    int numCols;
    std::vector<std::vector<ComplexNum>> entryData; //Outer vector is row, inner vector is column. For example, entryData[1][2] accesses the second row, third column
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

ComplexNum InnerProduct(std::vector<ComplexNum>* u, std::vector<ComplexNum>* v) {
    if (u->size() != v->size()) {
        throw std::invalid_argument("Invalid Dimensions");
    }

    ComplexNum result(0.0,0.0);
    for(int i = 0; i < u->size(); i++) {
        result += (*u)[i] * (*v)[i];
    }
    return result;
}

//improve once vector class is made will make all this stuff way smaller

Matrix GramSchmidt(Matrix* M) {
    Matrix result(M->numRows,M->numCols);

    for(int i = 1; i <= M->numCols; i++) {
        std::vector<ComplexNum> Vk = M->operator[](i);
        std::vector<ComplexNum> Uk = Vk;

        for(int j = 1; j < i; j++) {

            std::vector<ComplexNum> Un = result[j];
            ComplexNum VdotU = InnerProduct(&Vk,&(Un));
            ComplexNum UdotU = InnerProduct(&(Un),&(Un));

            //can make this read better with our own vector class functions for scalar multiplication
            ComplexNum quotient = (VdotU / UdotU);
            std::vector<ComplexNum> projection = Un;
            for(int k = 0; k < projection.size(); k++) {
                projection[k] = projection[k] * quotient;
            }

            for(int k = 0; k < projection.size(); k++) {
                Uk[k] = Uk[k] - projection[k];
            }
        }

        result.columnAssign(i,&Uk);
        std::cout << "here";
    }
    std::cout << "here";
    return result;
}

int main() {
    // Complex number test asserts
    ComplexNum z(3, 2);
    assert(z.getConjugate().complexPart == -2 );
    assert ((2 + z).realPart == 5);
    assert((z * -7).complexPart == -14);
    assert((2 * z).realPart == 6);
    assert((z * 12).realPart == 36);
    assert((z * 0).realPart == 0);
    assert(z.getMagnitude() == sqrt(13));
    ComplexNum w(2, 3);
    assert((z*w).complexPart == 13);
    assert((z*w).realPart == 0);



    // Matrix test asserts
    Matrix A(2, 2);
    A(1,2) = ComplexNum(1, 0);
    assert(A(1, 2) == ComplexNum(1, 0));
    Matrix B(2, 2);
    B(1, 2) = ComplexNum(3, 4);
    Matrix C = A + B;
    assert(C(1, 2) == ComplexNum(4, 4));

    Matrix D(2,2);
    D(1, 1) = 5;
    D(1,2) = 1;
    D(2,1) = 1;
    D(2, 2) = 6;
    auto G = GramSchmidt(&D);
    std::cout << "dnak";






};
#pragma clang diagnostic pop