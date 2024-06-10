#include "ComplexNum.h"


// Constructor
ComplexNum::ComplexNum(double a, double b) : realPart(a), imagPart(b) {}

// Equality operators
ComplexNum& ComplexNum::operator=(ComplexNum const& numToCopy) {
    if (this != &numToCopy) {
        this->realPart = numToCopy.realPart;
        this->imagPart = numToCopy.imagPart;
    }
    return *this;
}

ComplexNum& ComplexNum::operator=(double numToCopy) {
    this->realPart = numToCopy;
    this->imagPart = 0;
    return *this;
}

ComplexNum& ComplexNum::operator+=(ComplexNum const& numToAdd) {
    this->realPart += numToAdd.realPart;
    this->imagPart += numToAdd.imagPart;
    return *this;
}

bool ComplexNum::operator==(const ComplexNum& otherNum) const {
    return (this->realPart == otherNum.realPart && this->imagPart == otherNum.imagPart);
}
bool ComplexNum::operator==(double otherNum) const {
    return (this->realPart == otherNum && this->imagPart == 0);
}
bool ComplexNum::operator!=(const ComplexNum& otherNum) const {
    return (!(this->realPart == otherNum.realPart && this->imagPart == otherNum.imagPart));
}

bool ComplexNum::operator<(const ComplexNum& otherNum) const {
    if (this->realPart < otherNum.realPart) {
        return true;
    }
    if (this->realPart > otherNum.realPart) {
        return false;
    }
    return this->imagPart < otherNum.imagPart;
}

// Copy constructors

ComplexNum::ComplexNum(const ComplexNum& numToCopy) {
    this->realPart = numToCopy.realPart;
    this->imagPart = numToCopy.imagPart;
}
ComplexNum::ComplexNum(double numToCopy) {
    this->realPart = numToCopy;
    this->imagPart = 0;
}

// Arithmetic operators

ComplexNum ComplexNum::operator+(ComplexNum const& numToAdd) const {
    return {this->realPart + numToAdd.realPart, this->imagPart + numToAdd.imagPart};
}
ComplexNum ComplexNum::operator+(double numToAdd) const {
    return {this->realPart + numToAdd, this->imagPart};
}

ComplexNum ComplexNum::operator*(ComplexNum const& numToMul) const {
    double a = this->realPart;
    double b = this->imagPart;
    double c = numToMul.realPart;
    double d = numToMul.imagPart;

    return {(a * c) - (b * d), (a * d) + (b * c)};
}

ComplexNum ComplexNum::operator*(double numToMul) const {
    return {numToMul * this->realPart, numToMul * this->imagPart};
}

ComplexNum ComplexNum::operator-(ComplexNum const& numToSub) const {
    return {this->realPart - numToSub.realPart, this->imagPart - numToSub.imagPart};
}
ComplexNum ComplexNum::operator-(double numToSub) const {
    return {this->realPart - numToSub, this->imagPart};
}

ComplexNum ComplexNum::operator/(double numToDiv) const {
    return {this->realPart / numToDiv, this->imagPart / numToDiv};
}
ComplexNum ComplexNum::operator/(ComplexNum const& numToDiv) const {
    double a2SquaredPlusb2Squared = (numToDiv.realPart * numToDiv.realPart) + (numToDiv.imagPart * numToDiv.imagPart);
    double a = (((this->realPart * numToDiv.realPart) + (this->imagPart * numToDiv.imagPart)) / a2SquaredPlusb2Squared);
    double b = ((this->imagPart * numToDiv.realPart) - (this->realPart * numToDiv.imagPart)) / a2SquaredPlusb2Squared;
    return {a,b};
}

std::ostream& operator<<(std::ostream& outputStream, const ComplexNum& numberToPrint) {
    double rp = numberToPrint.getRealPart();
    double ip = numberToPrint.getImagPart();

    if (rp != 0) {
        outputStream << rp;
    }
    if (ip != 0) {
        if (ip > 0 && rp != 0) {
            outputStream << "+";
        }
        if (ip < 0) {
            outputStream << "-";
            ip = -ip; // Make imaginary part positive for output
        }
        outputStream << ip << 'i';
    }
    if (rp == 0 && ip == 0) {
        outputStream << "0";
    }
    return outputStream;
}

//Accessors
double ComplexNum::getRealPart() const {
    return this->realPart;
}
double ComplexNum::getImagPart() const {
    return this->imagPart;
}

//Mutators
void ComplexNum::setRealPart(double value) {
    this->realPart = value;
}
void ComplexNum::setImagPart(double value) {
    this->imagPart = value;
}

//Methods
double ComplexNum::getMagnitude() const {
    if (imagPart == 0) {
        return abs(realPart);
    }
    return sqrt((realPart * realPart) + (imagPart * imagPart));
}
ComplexNum ComplexNum::getConjugate() const {
    ComplexNum numToReturn;
    numToReturn.realPart = this->realPart;
    numToReturn.imagPart = this->imagPart * -1;
    return numToReturn;
}

double ComplexNum::sign() const {
    if(this->realPart == 0 || this->imagPart == 0) {
        return 1;
    }
    return ((*this) / magnitudeOfNumber(*this)).getRealPart();

}

// Friend functions for commutativity
ComplexNum operator+(double lhs, const ComplexNum& rhs) {
    return {lhs + rhs.getRealPart(), rhs.getImagPart()};
}
ComplexNum operator*(double lhs, const ComplexNum& rhs) {
    return {lhs * rhs.getRealPart(), lhs * rhs.getImagPart()};
}


//External funtions for easier use
ComplexNum complexConjugate(const ComplexNum& numToConjugate) {
    ComplexNum numToReturn;
    numToReturn.setRealPart(numToConjugate.getRealPart());
    numToReturn.setImagPart(numToConjugate.getImagPart() * -1);
    return numToReturn;
}
ComplexNum complexConjugate(double numToConjugate) {
    ComplexNum numToReturn;
    numToReturn.setRealPart(numToConjugate);
    numToReturn.setImagPart(0);
    return numToReturn;
}
double magnitudeOfNumber(const ComplexNum& numToMag) {
    return sqrt((numToMag.getRealPart() * numToMag.getRealPart()) + (numToMag.getImagPart() * numToMag.getImagPart()));
}
double magnitudeOfNumber(double numToMag) {
    return abs(numToMag);
}
double complexNumToDouble(const ComplexNum& num) {
    if (num.getImagPart() != 0) {
        throw std::invalid_argument("complexNumToDouble: Cannot convert complex number with non-zero imaginary part to double");
    }
    return num.getRealPart();
}