#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>


// Classes

class ComplexNum {
public:
    // Constructor
    ComplexNum(double a = 0, double b = 0) : realPart(a), complexPart(b) {};

    // Equality operators
    ComplexNum operator=(ComplexNum const& numToCopy) {
        if (this != &numToCopy) {
            this->realPart = numToCopy.realPart;
            this->complexPart = numToCopy.complexPart;
        }
        return *this;
    }
    ComplexNum operator=(double numToCopy) {
        this->realPart = numToCopy;
        this->complexPart = 0;
        return *this;
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
    // todo(?) division of two complex numbers

    // Method declarations
    double getMagnitude() {
        return sqrt((realPart * realPart) + (complexPart * complexPart));
    };
    double getRealPart() {
        return realPart;
    }
    double getComplexPart() {
        return complexPart;
    }
    ComplexNum getConjugate() {
        ComplexNum numToReturn;
        numToReturn.realPart = this->realPart;
        numToReturn.complexPart = this->complexPart * -1;
        return numToReturn;
    }

    // Internal values
    double realPart;
    double complexPart;
};










int main() {
    // Complex number test asserts
    ComplexNum z(3, 2);
    assert(z.getConjugate().complexPart == -2 );
    assert ((z + 2).realPart == 5);
    assert((z * -7).complexPart == -14);
    assert((z * 2).realPart == 6);
    assert((z * 12).realPart == 36);
    assert((z * 0).realPart == 0);
    assert(z.getMagnitude() == sqrt(13));
    ComplexNum w(2, 3);
    assert((z*w).complexPart == 13);
    assert((z*w).realPart == 0);

};