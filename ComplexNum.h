#pragma once
#include <math.h>
#include <stdexcept>
#include <iostream>

class ComplexNum {
private:
    //Internal values
    double realPart;
    double imagPart;
public:
    //Constructor
    ComplexNum (double a = 0, double b= 0);

    //Equality operators
    ComplexNum& operator=(ComplexNum const& numToCopy);
    ComplexNum& operator=(double numToCopy);
    ComplexNum& operator+=(ComplexNum const& numToAdd);
    bool operator==(const ComplexNum& otherNum) const;
    bool operator==(double otherNum) const;
    bool operator!=(const ComplexNum& otherNum) const;

    //Copy constructors
    ComplexNum(const ComplexNum& numToCopy);
    ComplexNum(double numToCopy);

    //Arithmetic operators
    ComplexNum operator+(ComplexNum const& numToAdd) const;
    ComplexNum operator+ (double numToAdd) const;
    ComplexNum operator*(ComplexNum const& numToMul) const;
    ComplexNum operator*(double numToMul) const;
    ComplexNum operator-(ComplexNum const& numToSub) const;
    ComplexNum operator-(double numToSub) const;
    ComplexNum operator/(double numToDiv) const;
    ComplexNum operator/(ComplexNum const& numToDiv) const;

    //Accessors
    double getRealPart() const;
    double getImagPart() const;

    //Mutators
    void setRealPart(double value);
    void setImagPart(double value);

    //Methods
    double getMagnitude() const;
    ComplexNum getConjugate() const;

    friend std::ostream& operator<<(std::ostream& outputStream, const ComplexNum& numberToPrint);

};

//commutivity
ComplexNum operator+(double lhs, const ComplexNum& rhs);
ComplexNum operator*(double lhs, const ComplexNum& rhs);

// External complex conjugate methods added for easier use
ComplexNum complexConjugate(const ComplexNum& numToConjugate);
ComplexNum complexConjugate(double numToConjugate);
double magnitudeOfNumber(const ComplexNum& numToMag);
double magnitudeOfNumber(double numToMag);
double complexNumToDouble(const ComplexNum& num);


