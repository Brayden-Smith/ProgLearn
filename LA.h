#pragma once
#import "Matrix.h"

Matrix matMul(Matrix const& lhs, Matrix const& rhs);
double realTraceOfMatrix(Matrix const& matrixToTrace);
Matrix conjTranspose(Matrix const& matrixToTranspose);
double frobeniusNorm(Matrix const& matrixToNorm);
ComplexNum InnerProduct(Matrix& u, Matrix& v);
void normalizeVectorsInMatrix(Matrix* pointerToMatrix);
Matrix GramSchmidt(Matrix const& M);


