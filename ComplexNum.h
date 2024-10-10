#pragma once
#include <math.h>
#include <stdexcept>
#include <iostream>

class ComplexNum {
private:
    //Cartesian coordinates
    double realPart;
    double imagPart;

    //Polar coordinates
    double magnitude;
    double theta;
    void calculateMagnitude();
    void calculateTheta();
    // !WARNING! YOU ARE EXPECTED TO CALCULATE POLAR COORDINATES IN ALL FUNCTIONS THAT NEED THEM !WARNING!
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
    bool operator<(const ComplexNum& otherNum) const; // DO NOT USE; THIS ONLY EXISTS TO MAKE STD::SET WORK

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
    ComplexNum operator^(ComplexNum const& degree);

    //Accessors
    double getRealPart() const;
    double getImagPart() const;
    double getMagnitude();
    double getTheta() const;

    //Mutators
    void setRealPart(double value);
    void setImagPart(double value);

    //Methods

    ComplexNum getConjugate() const;
    double sign() const;

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
ComplexNum sqrt(ComplexNum numToSquareroot);
double getRealSign(ComplexNum num);

ComplexNum power(ComplexNum num, int power);


