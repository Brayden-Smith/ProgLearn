#pragma once
#include "Matrix.h"
#include <set>
#include <unordered_set>
#include <opencv2/opencv.hpp>
#include "data.h"


void thresholdStabilize(Matrix* matrixToStabilize);

Matrix matMul(Matrix* lhs, Matrix* rhs);
Matrix hadamardProduct(Matrix* lhs, Matrix* rhs);

double realTraceOfMatrix(Matrix* matrixToTrace);
Matrix conjTranspose (Matrix* matrixToTranspose);
Matrix transpose (Matrix* matrixToTranspose);
ComplexNum frobeniusNorm(Matrix* matrixToNorm);

ComplexNum innerProduct(Matrix& u, Matrix& v);
ComplexNum mean(Matrix& u);

void normalizeVectorsInMatrix(Matrix* pointerToMatrix);
Matrix GramSchmidt(Matrix const& M);
Matrix createAugmentedMatrix(Matrix* matrixA, Matrix* matrixB);
void swapRowsInMatrix(Matrix* matrixToSwap, int posRowOne, int posRowTwo);
Matrix gaussianElimination(Matrix* matrix, Matrix* vector);
Matrix inverseMatrix(Matrix* matrixToInvert);

//needed for QR
Matrix minorMatrix(Matrix M,int iDel, int jDel);
double VectorNorm(Matrix* pointerToMatrix);
Matrix unitVector(int k, int dim);
Matrix identityMatrix(int dim);
Matrix frontFillVec(Matrix M, int dim, ComplexNum const& Fill);
Matrix tridiagonalizeMatrix(Matrix& matrixToTri); // For our beloved symmetric AA* matrices
Matrix householderReflection(Matrix* x);
std::vector<Matrix> QRDecomp(Matrix const& A);

bool isUpperTriangular(Matrix* mat);
std::vector<ComplexNum> eigenvalues(Matrix* matrix);
void orthogonalize(Matrix* Q);
std::vector<Matrix> eigenvectors(Matrix* matrix, std::vector<ComplexNum>& correspondingEigenValues);
std::vector<Matrix> singularValueDecomp(Matrix* matrix);


ComplexNum expectedValue(Matrix* Z);
ComplexNum covariance(Matrix* Z, Matrix* W);
Matrix covarianceMatrix(Matrix* M);

Matrix littleCovariance(Matrix& matrix);

std::vector<ComplexNum> CVEigenValues(Matrix& matrix); // Only when real eigenvalues are expected (symmetric matrices)
std::vector<Matrix> CVEigenvectors(Matrix* matrix, std::vector<ComplexNum>& correspondingEigenValues); // Only when real eigenvalues are expected (symmetric matrices)

double smallestNumberInRealMatrix(Matrix& matrix);
double largestNumberInRealMatrix(Matrix& matrix);
Matrix normalizeMatrix(Matrix& matrix, double maxVal, double minVal);