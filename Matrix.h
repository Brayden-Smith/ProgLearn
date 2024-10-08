#pragma once
#include "ComplexNum.h"
#include <vector>
#include <iostream>

class Matrix  {
private:
    int numRows;
    int numCols;
    std::vector<std::vector<ComplexNum>> entryData;
public:
    //Constructors
    Matrix (int nr, int nc);
    Matrix(const Matrix& matToCopy);
    Matrix();

    //Arithmetic operators
    Matrix operator + (Matrix const& matrixToAdd) const;

    Matrix operator *(ComplexNum const& scalar);
    Matrix operator *(double numToMul);
    friend Matrix operator*(ComplexNum const& scalar, Matrix const& matrix);
    friend Matrix operator*(double numToMul, Matrix const& matrix);

    Matrix operator /(ComplexNum const& numToDiv);

    Matrix operator -(Matrix const& matrixToSub) const;
    Matrix operator -(ComplexNum const& numToSub) const;


    //Equality operators
    Matrix& operator =(Matrix const& matrixToCopy);
    bool operator ==(const Matrix& otherMatrix) const;

    //Entry accessor operator
    ComplexNum& operator ()(int i, int j);

    //Column accessor operator
    Matrix operator[](int n) const;

    //Column mutators
    void columnAssign(int n, Matrix* M);

    //Row access operator
    Matrix operator()(int n) const;

    //Row mutator
    Matrix rowAssign(int n, Matrix* M);

    //Accessor methods
    int getNumRows() const;
    int getNumCols() const;

    friend std::ostream& operator<<(std::ostream& outputStream, const Matrix& matrixToPrint);

    Matrix conjugate() const;



};

bool hasNaN(Matrix& mat);
