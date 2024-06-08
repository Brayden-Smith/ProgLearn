#pragma once
#include "Matrix.h"

void thresholdStabilize(Matrix* matrixToStabilize);
Matrix matMul(Matrix* lhs, Matrix* rhs);
double realTraceOfMatrix(Matrix* matrixToTrace);
Matrix conjTranspose (Matrix* matrixToTranspose);
Matrix transpose (Matrix* matrixToTranspose);
double frobeniusNorm(Matrix* matrixToNorm);
ComplexNum InnerProduct(Matrix& u, Matrix& v);
void normalizeVectorsInMatrix(Matrix* pointerToMatrix);
Matrix GramSchmidt(Matrix const& M);
Matrix createAugmentedMatrix(Matrix* matrixA, Matrix* matrixB);
void swapRowsInMatrix(Matrix* matrixToSwap, int posRowOne, int posRowTwo);
Matrix gaussianElimination(Matrix* matrix, Matrix* vector);
Matrix inverseMatrix(Matrix* matrixToInvert);