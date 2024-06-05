#pragma once
#include "Matrix.h"

Matrix matMul(Matrix* lhs, Matrix* rhs);
double realTraceOfMatrix(Matrix* matrixToTrace);
Matrix conjTranspose(Matrix* matrixToTranspose);
double frobeniusNorm(Matrix* matrixToNorm);
ComplexNum InnerProduct(Matrix& u, Matrix& v);
void normalizeVectorsInMatrix(Matrix* pointerToMatrix);
Matrix GramSchmidt(Matrix const& M);